import re

def parse_expression(expr, num_vars):
    coeffs = [0.0] * num_vars
    terms = re.findall(r'([+-]?\s*\d*\.?\d*/?\d*\*?x\d+)', expr.replace(" ", ""))
    for term in terms:
        match = re.match(r'([+-]?\d*\.?\d*/?\d*)\*?x(\d+)', term)
        if match:
            coef_str, var_idx = match.groups()
            coef = 1.0 if coef_str in ["", "+"] else -1.0 if coef_str == "-" else eval(coef_str)
            coeffs[int(var_idx)-1] = coef
    return coeffs

def ler_sistema_arquivo(nome_arquivo):
    with open(nome_arquivo, 'r') as f:
        linhas = [linha.strip() for linha in f if linha.strip()]

    func_obj = linhas[0].lower()
    is_min = "min" in func_obj
    num_vars = max(map(int, re.findall(r'x(\d+)', func_obj)))

    C = parse_expression(func_obj.split('=')[1], num_vars)
    if not is_min:
        C = [-c for c in C]

    A, b, tipos = [], [], []
    for restr in linhas[1:]:
        for tipo in ["<=", ">=", "="]:
            if tipo in restr:
                lhs, rhs = restr.split(tipo)
                A.append(parse_expression(lhs, num_vars))
                b.append(float(eval(rhs.strip())))
                tipos.append(tipo)
                break

    for i, tipo in enumerate(tipos):
        for j in range(len(A)):
            A[j].append(1 if j == i and tipo == "<=" else -1 if j == i and tipo == ">=" else 0)

    while len(C) < len(A[0]):
        C.append(0.0)

    return A, b, C, is_min

def determinante_laplace(matriz):
    n = len(matriz)
    if n == 1: return matriz[0][0]
    if n == 2: return matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0]
    return sum(((-1)**j) * matriz[0][j] * determinante_laplace([linha[:j] + linha[j+1:] for linha in matriz[1:]]) for j in range(n))

def inversa_matriz(matriz):
    det = determinante_laplace(matriz)
    if abs(det) < 1e-8:
        raise ValueError("Matriz não invertível (det = 0).")
    cofatores = [[((-1)**(i+j)) * determinante_laplace([linha[:j] + linha[j+1:] for k, linha in enumerate(matriz) if k != i]) for j in range(len(matriz))] for i in range(len(matriz))]
    adjunta = [list(col) for col in zip(*cofatores)]
    return [[elem / det for elem in linha] for linha in adjunta]

def precisa_de_fase_I(A, b):
    # Verifica se há alguma linha da matriz A sem 1 como entrada básica clara
    num_restricoes = len(A)
    num_variaveis = len(A[0])
    for col in range(num_variaveis):
        coluna = [A[row][col] for row in range(num_restricoes)]
        if coluna.count(1) == 1 and coluna.count(0) == num_restricoes - 1:
            continue
        else:
            return True
    return any(bi < 0 for bi in b)

def fase_I(A, b):
    A_aux = [linha[:] for linha in A]  # cópia manual (lista de listas)
    artificiais = []
    for i in range(len(A)):
        for j in range(len(A)):
            A_aux[j].append(1.0 if j == i else 0.0)
        artificiais.append(len(A_aux[0]) - 1)
    C_fase1 = [0.0] * (len(A_aux[0]) - len(artificiais)) + [1.0] * len(artificiais) + [0.0]
    base = artificiais[:]
    tabela = [linha + [b[i]] for i, linha in enumerate(A_aux)]

    while True:
        Z = [sum(C_fase1[base[j]] * tabela[j][i] for j in range(len(base))) - C_fase1[i] for i in range(len(C_fase1))]
        if all(z <= 1e-8 for z in Z[:-1]): break
        col_pivo = Z.index(max(Z[:-1]))
        razoes = [tabela[i][-1] / tabela[i][col_pivo] if tabela[i][col_pivo] > 1e-8 else float('inf') for i in range(len(A))]
        linha_pivo = razoes.index(min(razoes))
        if all(r == float('inf') for r in razoes):
            raise Exception("Solução ilimitada.")
        base[linha_pivo] = col_pivo
        pivo = tabela[linha_pivo][col_pivo]
        tabela[linha_pivo] = [x / pivo for x in tabela[linha_pivo]]
        for i in range(len(tabela)):
            if i != linha_pivo:
                fator = tabela[i][col_pivo]
                tabela[i] = [tabela[i][j] - fator * tabela[linha_pivo][j] for j in range(len(C_fase1))]

    # Verificar se há variáveis artificiais ainda na base com valor diferente de zero
    for i, var in enumerate(base):
        if var in artificiais and abs(tabela[i][-1]) > 1e-8:
            raise Exception("Problema inviável. Valor ótimo da Fase I > 0 (variável artificial ativa).")


    for idx in sorted(artificiais, reverse=True):
        for linha in tabela:
            del linha[idx]
    base = [v for v in base if v not in artificiais]

    return tabela, base

def fase_II(tabela, base, C_original):
    num_vars = len(C_original)
    C = C_original + [0.0] * (len(tabela[0]) - num_vars - 1) + [0.0]

    while True:
        Z = [sum(C[base[j]] * tabela[j][i] for j in range(len(base))) - C[i] for i in range(len(C))]
        if all(z <= 1e-8 for z in Z[:-1]): break
        col_pivo = Z.index(max(Z[:-1]))
        razoes = [tabela[i][-1] / tabela[i][col_pivo] if tabela[i][col_pivo] > 1e-8 else float('inf') for i in range(len(tabela))]
        linha_pivo = razoes.index(min(razoes))
        if razoes[linha_pivo] == float('inf'):
            raise Exception("Solução ilimitada.")
        base[linha_pivo] = col_pivo
        pivo = tabela[linha_pivo][col_pivo]
        tabela[linha_pivo] = [x / pivo for x in tabela[linha_pivo]]
        for i in range(len(tabela)):
            if i != linha_pivo:
                fator = tabela[i][col_pivo]
                tabela[i] = [tabela[i][j] - fator * tabela[linha_pivo][j] for j in range(len(C))]

    solucao = [0.0] * num_vars
    for i, var in enumerate(base):
        if var < num_vars:
            solucao[var] = tabela[i][-1]
    valor_otimo = sum(C_original[i] * solucao[i] for i in range(num_vars))
    return solucao, valor_otimo

def simplex(A, b, C):
    print("[INFO] Iniciando Simplex...")
    if precisa_de_fase_I(A, b):
        print("[INFO] Executando Fase I...")
        tabela, base = fase_I(A, b)
    else:
        print("[INFO] Executando Fase II diretamente...")
        tabela = [linha[:] + [b[i]] for i, linha in enumerate(A)]
        base = []
        for col in range(len(A[0])):
            coluna = [A[row][col] for row in range(len(A))]
            if coluna.count(1) == 1 and coluna.count(0) == len(A) - 1:
                base.append(col)
                if len(base) == len(A):
                    break
    return fase_II(tabela, base, C)

if __name__ == "__main__":
    A, b, C, is_min = ler_sistema_arquivo("entrada.txt")
    print("\n[INFO] Sistema lido:")
    print("Matriz A:")
    for linha in A:
        print(linha)
    print("\nVetor b:", b)
    print("\nVetor C (função objetivo):", C)

    try:
        solucao, valor_otimo = simplex(A, b, C)
        if not is_min:
            valor_otimo *= -1
        print("\n[RESULTADO FINAL]")
        print("Solução ótima:", solucao)
        print("Valor ótimo da função objetivo:", valor_otimo)
    except Exception as e:
        print("[ERRO]", e)
