def reversal_distance_manual(p1, p2):
    """
    Вычисляет расстояние обращения (reversal distance) между двумя перестановками
    с использованием алгоритма поиска в ширину (BFS).
    
    Параметры:
        p1 (list): Исходная перестановка.
        p2 (list): Целевая перестановка.
    
    Возвращает:
        int: Минимальное количество обращений для преобразования p1 в p2.
    """
    # Нормализуем перестановки: преобразуем p2 в [1,2,3,...], а p1 соответственно
    normalized_index = {value: i + 1 for i, value in enumerate(p2)}
    normalized_start = [normalized_index[value] for value in p1]
    target = list(range(1, len(p1) + 1))
    
    if normalized_start == target:
        return 0
    
    queue = [(normalized_start, 0)]
    visited = {tuple(normalized_start)}
    
    while queue:
        current_perm, dist = queue.pop(0)
        
        # Генерируем все возможные обращения для текущей перестановки
        for i in range(len(current_perm)):
            for j in range(i + 1, len(current_perm)):
                new_perm = current_perm[:i] + list(reversed(current_perm[i:j+1])) + current_perm[j+1:]
                new_perm_tuple = tuple(new_perm)
                
                if new_perm_tuple not in visited:
                    if new_perm == target:
                        return dist + 1
                    visited.add(new_perm_tuple)
                    queue.append((new_perm, dist + 1))
    
    return -1  # На случай, если решение не найдено (для перестановок одинаковой длины должно быть найдено всегда)

# --- Функции для чтения данных и вывода результата ---
def read_permutations_from_string(data_string):
    """
    Читает пары перестановок из строки в формате, аналогичном файлу Rosalind.
    
    Параметры:
        data_string (str): Входные данные в виде строки.
    
    Возвращает:
        list: Список пар перестановок (каждая пара - кортеж из двух списков).
    """
    lines = data_string.strip().split('\n')
    pairs = []
    i = 0
    while i < len(lines):
        if lines[i].strip():
            # Читаем две последовательные строки как пару перестановок
            perm1 = list(map(int, lines[i].split()))
            perm2 = list(map(int, lines[i + 1].split()))
            pairs.append((perm1, perm2))
            i += 2
        else:
            i += 1
    return pairs

def main():
    """
    Основная функция для чтения входных данных из файла, вычисления расстояний
    и вывода результатов.
    """
    import sys
    
    # Чтение данных из файла (стандартный ввод или файл 'rosalind_rear.txt')
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'r') as f:
            data_string = f.read()
    else:
        try:
            with open('Dataset.txt', 'r') as f:
                data_string = f.read()
        except FileNotFoundError:
            # Если файл не найден, читаем со стандартного ввода
            data_string = sys.stdin.read()
    
    # Парсинг пар перестановок
    permutation_pairs = read_permutations_from_string(data_string)
    
    # Вычисление расстояния обращения для каждой пары
    distances = []
    for perm1, perm2 in permutation_pairs:
        distances.append(str(reversal_distance_manual(perm1, perm2)))
    
    # Вывод результатов через пробел
    print(' '.join(distances))

if __name__ == "__main__":
    main()