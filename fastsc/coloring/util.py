"""
util.py -  Utility functions
"""

"""
Input: A dictionary with (key, value) = (edge, int color)
Return: A dictionary with (key, value) = (edge, int color), such that if a<b
then more nodes are colored with a than b
"""
def relabel_coloring(int_coloring):
    num_int = len(set(int_coloring.values()))
    arr = [] # [(color, popularity) for each color]
    for i in range(num_int):
        num = len([k for k,v in int_coloring.items() if v == i])
        arr.append((i,num))
    arr.sort(key=lambda x:x[1]) # sort arr by popularity
    new_coloring = {}
    for i in range(num_int):
        old_color = arr[num_int-i-1][0] # the new color is i
        # collect the list of nodes with color = old_color
        temp = [k for k,v in int_coloring.items() if v == old_color]
        # color them with the new color
        for k in temp:
            new_coloring[k] = i
    return new_coloring
