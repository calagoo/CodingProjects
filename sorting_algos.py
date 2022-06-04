from array import array
import os
import sys
import random
import matplotlib.pyplot as plt

def rand_data(array_length):
    random_array = []
    while len(random_array) < array_length:
        random_int = random.randint(1,array_length)
        if random_int in random_array:
            continue
        else:
            random_array.append(random_int)
    return random_array

def bubble_sort(array):
    for j in range(0,len(array)-1):
        if array[j] > array[j+1]:
            array[j+1], array[j] = array[j], array[j+1]
    return array

def insertion_sort(array,speed):
    for i in range(1,len(array)):
        key = array[i]
        j = i - 1
        while j >= 0 and array[j] > key:
            array[j+1], array[j] = array[j], array[j+1]
            j = j - 1
        array[j+1] = key
        plt.clf()
        plt.bar(range(len(array)),array)
        plt.title("Insertion")
        plt.draw()
        plt.pause(1/speed)

def bar_plot(array,sort_type,speed):
    plt.figure()
    for _ in range(len(array)):
        if sort_type == 'bubble':
            arr = bubble_sort(array)
            plt.clf()
            plt.bar(range(len(arr)),arr)
            plt.title("Bubble")
            plt.draw()
            plt.pause(1/speed)
    if sort_type == 'insert':
        arr = insertion_sort(array,speed)

def main():
    arrlen = 50
    speed = 100
    arr = rand_data(array_length=arrlen)
    bar_plot(array=arr,sort_type='bubble',speed=speed)
    plt.show()

if __name__ == "__main__":
    main()