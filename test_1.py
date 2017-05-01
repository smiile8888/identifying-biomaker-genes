items_list = [24, 150, 79, 88, 345, 3, 50]
credit = 200

for i in range(len(items_list[:-1])):
    first = i
    last_index = len(items_list) - 1
    for j in range(i+1, last_index):
        if items_list[i] + items_list[j+1] == credit:
            second = j+1
            print first + 1
            print second + 1
            exit()
