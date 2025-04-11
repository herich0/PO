import re

def parse_expression(expr, num_vars):
    coeffs = [0.0] * num_vars
    terms = re.findall(r'([+-]?\s*\d*\.?\d*/?\d*\*?x\d+)', expr.replace(" ", ""))
    for term in terms:
        term = term.replace(" ", "")
        match = re.match(r'([+-]?\d*\.?\d*/?\d*)\*?x(\d+)', term)
        if match:
            coef_str, var_idx = match.groups()
            if coef_str in ["", "+"]:
                coef = 1.0
            elif coef_str == "-":
                coef = -1.0
            elif '/' in coef_str:
                numerador, denominador = coef_str.split('/')
                coef = float(numerador) / float(denominador)
            else:
                coef = float(coef_str)
            coeffs[int(var_idx)-1] = coef
    return coeffs


def ler_sistema_arquivo(nome_arquivo):
    with open(nome_arquivo, 'r') as f:
        linhas = [linha.strip() for linha in f if linha.strip()]

    func_obj = linhas[0].lower()
    is_min = "min" in func_obj
    match_vars = re.findall(r'x(\d+)', func_obj)
    num_vars = max(map(int, match_vars))

    C = parse_expression(func_obj.split('=')[1], num_vars)

    # Converte min em max (negando os coeficientes)
    if is_min:
        C = [-c for c in C]

    A = []
    b = []

    folga_index = 0  # número da folga atual
    tipos_restricoes = []  # guarda tipo: <=, >=, =

    for restr in linhas[1:]:
        if "<=" in restr:
            lhs, rhs = restr.split("<=")
            tipo = "<="
        elif ">=" in restr:
            lhs, rhs = restr.split(">=")
            tipo = ">="
        elif "=" in restr:
            lhs, rhs = restr.split("=")
            tipo = "="
        else:
            continue

        linha_A = parse_expression(lhs, num_vars)
        tipos_restricoes.append(tipo)
        b.append(int(rhs.strip()))
        A.append(linha_A)

    num_restricoes = len(A)

    folga_colunas = 0
    for i, tipo in enumerate(tipos_restricoes):
        if tipo in ["<=", ">="]:
            for j in range(num_restricoes):
                if j == i:
                    A[j].append(1 if tipo == "<=" else -1)
                else:
                    A[j].append(0)
            folga_colunas += 1

    # Ajustar C para ter o mesmo número de colunas que A
    while len(C) < len(A[0]):
        C.append(0.0)


    # Ajustar C para ter mesmo número de colunas que A
    while len(C) < len(A[0]):
        C.append(0)

    return A, b, C

def extrair_matriz_b(A):
    m = len(A)  # número de restrições (linhas)
    B = [linha[:m] for linha in A]
    return B

def determinante_laplace(matriz):
    n = len(matriz)

    # Caso base: 1x1
    if n == 1:
        return matriz[0][0]
    
    # Caso base: 2x2
    if n == 2:
        return matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0]
    
    # Expansão de Laplace na primeira linha
    det = 0
    for j in range(n):
        submatriz = [linha[:j] + linha[j+1:] for linha in matriz[1:]]
        cofator = ((-1) ** j) * matriz[0][j] * determinante_laplace(submatriz)
        det += cofator
    
    return det

def matriz_cofatores(matriz):
    n = len(matriz)
    cofatores = []

    for i in range(n):
        linha_cof = []
        for j in range(n):
            # Submatriz sem linha i e coluna j
            submatriz = [linha[:j] + linha[j+1:] for k, linha in enumerate(matriz) if k != i]
            cof = ((-1) ** (i + j)) * determinante_laplace(submatriz)
            linha_cof.append(cof)
        cofatores.append(linha_cof)

    return cofatores

def transpor_matriz(matriz):
    return [list(col) for col in zip(*matriz)]

def inversa_matriz(matriz):
    det = determinante_laplace(matriz)
    if det == 0:
        raise ValueError("A matriz não é invertível (determinante zero).")
    
    cofatores = matriz_cofatores(matriz)
    adjunta = transpor_matriz(cofatores)
    inversa = [[elem / det for elem in linha] for linha in adjunta]
    return inversa

# Teste
A, b, C = ler_sistema_arquivo("entrada.txt")
B = extrair_matriz_b(A)
det_B = determinante_laplace(B)

print("Matriz A:")
for linha in A:
    print(linha)

print("\nVetor b:")
print(b)

print("\nVetor C:")
print(C)

print("\nMatriz B:")
for linha in B:
    print(linha)

print("\nDeterminante de B:")
print(det_B)

try:
    inv_B = inversa_matriz(B)

    print("\nMatriz inversa de B:")
    for linha in inv_B:
        print(linha)
except ValueError as e:
    print(f"\nErro: {e}")
