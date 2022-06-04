import os
from statistics import mean
import sys
import random
import time
from sys import platform
from collections import OrderedDict

# Clear the system
def clear_system():
    if platform == "linux" or platform == "linux2":
        os.system("clear")
    elif platform == "win32":
        os.system("cls")

def wordle_wordlist(guess_words,answer_words):
    word_list = []
    with open(guess_words,'r') as file:
        lines = file.readlines()
        for word in lines:
            word_list.append(word.replace('\n',''))
    with open(answer_words,'r') as file:
        lines = file.readlines()
        for word in lines:
            word_list.append(word.replace('\n',''))
    return word_list

def random_word(word_list):
    wordle_answer = random.choice(word_list)
    return wordle_answer

def wordle_solver_pure_random_guess(word_list,wordle_answer):
    total_attempts = 6
    
    for attempt in range(total_attempts):
        print('Attempt =',attempt+1)
        wordle_guess = random.choice(word_list)
        print(wordle_guess)
        if wordle_guess == wordle_answer:
            return print('Congratulations! Somehow you guessed the correct answer!')
    return print("Aw drats, you didn't get it. Who woulda thought!")

def wordle_solver_semi_random_guess(first_guess,word_list,wordle_answer):

    wordle_answer_split = list(wordle_answer)
    wordle_guess_split_list = []
    matching_letters = []
    matching_letters_yellow = []
    unmatching_letters = []
    total_attempts = 6
    attempt = 0
    while True:
        # wordle_guess = random.choice(word_list)
        if attempt == 0:
            wordle_guess_split = first_guess
        else:
            wordle_guess_split = list(random.choice(word_list))
        # print(''.join(wordle_guess_split))
        # time.sleep(.1)
        ## Checking if matching letters are in the guessed word

        def checking_matched_letters(matching_letters,wordle_guess_split):
            if not matching_letters == []:
                for matched_letter in matching_letters:
                    if matched_letter[0] == wordle_guess_split[matched_letter[1]]:
                        all_matched = True
                    else:
                        all_matched = False
                        return '', False
                if all_matched:
                    final_guess = wordle_guess_split
                    return final_guess, True
                else:
                    return '', False

        if not unmatching_letters == []:
            unmatched_letters = any(elem in unmatching_letters for elem in wordle_guess_split)
            if unmatched_letters == True:
                continue

        if not matching_letters_yellow == []:
            matched_letters_yellow = all(elem in wordle_guess_split for elem in matching_letters_yellow)
            if matched_letters_yellow == False:
                continue

        if not matching_letters == []:
            final_guess_and_bool = checking_matched_letters(matching_letters,wordle_guess_split)
            if final_guess_and_bool[1] == False:
                continue
            else:
                wordle_guess_split = final_guess_and_bool[0]
        if wordle_guess_split in wordle_guess_split_list:
            continue
        ## Checking if any letters match
        letter_index = 0
        for letter in wordle_guess_split:
            if letter in wordle_answer_split:
                # This is a yellow letter, all green letters are also yellow letters
                matching_letters_yellow.append(letter)
                if letter == wordle_answer_split[letter_index]:
                    # This would mean that the letters at a certain index match (aka a green letter in wordle)
                    # print("Letter at index", letter_index, "matches")
                    matching_letters.append([letter,letter_index])
            elif letter not in wordle_answer_split:
                unmatching_letters.append(letter)
            letter_index += 1

        wordle_guess_split_list.append(wordle_guess_split)
        ## Removing duplicates from matching_letters using a reserved list and list comprehension
        reserve_list1 = []
        reserve_list2 = []
        reserve_list3 = []
        [reserve_list1.append(x) for x in matching_letters if x not in reserve_list1]
        [reserve_list2.append(x) for x in matching_letters_yellow if x not in reserve_list2]
        [reserve_list3.append(x) for x in unmatching_letters if x not in reserve_list3]
        matching_letters = reserve_list1
        matching_letters_yellow = reserve_list2
        unmatching_letters = reserve_list3

        # matching_letters
        # print(matching_letters)
        # print(matching_letters_yellow)
        # print(unmatching_letters)

        ## checking win and lose conditions
        if wordle_guess_split_list[attempt] == wordle_answer_split:
            return 'Congratulations! Somehow you guessed the correct answer!', True, [''.join(x) for x in wordle_guess_split_list],attempt+1

        attempt += 1
        # print('Attempt =',attempt,'\n')
        if attempt >= 6:
            return "Aw drats, you didn't get it. Who woulda thought!", False, [''.join(x) for x in wordle_guess_split_list],6


if __name__ == "__main__":


    ## Change path to file location
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    # File names of word lists
    words_file_path = "common_words"
    target_words_file_path = "target_words"

    # Functions start
    clear_system()
    word_list = wordle_wordlist(words_file_path,target_words_file_path)
    wordle_answer = random_word(word_list)

    print(f"{wordle_answer=}")
    # wordle_solver_pure_random_guess(word_list,wordle_answer)
    first_guess = 'crane'
    results = wordle_solver_semi_random_guess(first_guess,word_list,wordle_answer)
    print(results[0])
    print('Words used: ' + str(results[2]))

    # # Uncomment to run and view the average amount of games taken over 100 Trials
    # games = []
    # attempts = []
    # attempt_sum = 0
    # first_guess = 'crane'
    # for Trials in range(100):
    #     game = 0
    #     while True:
    #         game += 1
    #         RESULTS = wordle_solver_semi_random_guess(first_guess,word_list,wordle_answer)

    #         # print(RESULTS[0])
    #         attempt_sum += RESULTS[3]
    #         if RESULTS[1]:
    #             # print(RESULTS[0])
    #             # print('It took',game,'games to reach the correct answer')
    #             games.append(game)
    #             attempts.append(attempt_sum)#RESULTS[3])
    #             attempt_sum = 0
    #             # print(RESULTS[3])
    #             break
    # print('On average, it takes',mean(games),'games to find the answer.\nThe average attempts are',mean(attempts))