import os
from string import punctuation
from sys import platform
from clear_system import clear_system
import time
from collections import Counter,OrderedDict
from string import punctuation
import re

# Clears terminal
clear_system()


def readSixLetter(txt_file,word_length):
    six_word = []
    with open(txt_file,"r",encoding="utf-8") as f:

        # reads whole document and splits by line
        lines = f.read().splitlines()

        # for loop that iterates line by line
        for line in lines:
            # replaces the punctuation with a space
            res = re.sub(r'[^\w\s]', ' ', line)
            # splits string by space
            words = res.split(' ')
            # for loop that iterates by each word in each line
            for word in words:
                # checks length of word and adds correct words to list
                if len(word.lower()) >= word_length:
                    six_word.append(word.lower())
    
    # finds and sorts frequency of words
    counts = Counter(six_word)
    new_list = sorted(six_word, key=lambda x: -counts[x])
    freq = Counter(new_list)

    return six_word, freq

if __name__ == "__main__":
    start_time = time.time()
    txt_file = "E:/Python/non_work/Ebook_Reader/don_quixote"
    word_length = 6
    six_word,freq = readSixLetter(txt_file,word_length)
    print("There are",len(six_word),"total words over 6 characters in length")
    indx = 0
    for key,val in list(freq.items())[0:10]:
        indx += 1
        print("{} - There are {} of the word '{}'".format(indx,val,key))
    print("It took {} seconds to run this code.".format(round(time.time()-start_time,4)))