"""
Script containing functionality related to DNA in Genetics like crossover & mutation.
Contains: 'do crossover' & 'do mutation' function
Author: Yasim Ahmad(yaaximus)
Email: yasim.ahmed63@yahoo.com
"""

from config import Config

import numpy as np
import random


def dna(chr_pop_fitness, ranked_population, chr_best_fitness_index, last_pop):
    """
    This function encapsulates functionality related to dna like crossover
    and mutation.
    Parameters
    ----------
    chr_pop_fitness : [numpy.ndarray]
        [Contains fitness values of chromosome population]
        [Chứa các giá trị thể chất của quần thể nhiễm sắc thể]
    ranked_population : [numpy.ndarray]
        [Contains numpy array of ranked chromosome population]
        [Chứa mảng numpy của quần thể nhiễm sắc thể được xếp hạng]
    chr_best_fitness_index : [list]
        [Contains list of best fitness indices in chromosome population]
        [Chứa danh sách các chỉ số thể dục tốt nhất trong quần thể nhiễm sắc thể]
    last_pop : [numpy.ndarray]
        [Contains numpy array of last population]
        [Chứa mảng của dân số cuối cùng]
    Returns
    -------
    [numpy.ndarray]
        [numpy array of chromosome with have gone through random crossover and mutation]
    """

    chromo_crossover_pop = _do_crossover(
        ranked_pop=ranked_population, chr_best_fit_indx=chr_best_fitness_index,
        pop=last_pop) #quá trình giao phối ADN

    chromo_crossover_mutated_pop = _do_mutation(pop=chromo_crossover_pop)#giao phối sau khi đột biến

    return chromo_crossover_mutated_pop


def _do_mutation(pop):#tạo đột biến
    """
    This function is responsible for handling mutation in population of chromosomes.
    Parameters
    ----------
    pop : [numpy.ndarray]
        [numpy array of chromosome population which will undergo mutation]
        [Mảng chứa các nhiễm sắc thể sẽ trải qua đột biến]
    Returns
    -------
    [numpy.ndarray]
        [numpy array of chromosome population undergone mutation]
        [Mảng chứa NST bị đột biến]
    """

    mutated_pop = np.array(pop, copy=True)#mảng NST đột biến

    itr = 3
    while itr < Config.pop_max:
        for k in range(Config.chr_len):
            c = random.random()
            if c < Config.mutation_rate and k != 0:
                mutated_pop[itr, k] = random.randint(1, Config.npts - 2)#Nếu tỉ lệ đột biến bé hơn 0.01 thì thay thế NST cũ bằng
                #NST mới trừ 2 điểm đầu và cuối
            else:
                pass
        itr += 1
    return mutated_pop


def _do_crossover(ranked_pop, chr_best_fit_indx, pop):#giao phối
    """
    This function is responsible for handling crossover in population of chromosomes.
    Parameters
    ----------
    ranked_pop : [numpy.ndarray]
        [numpy array of chromosome population which will undergo crossover]
        [Mảng chứa các pop sẽ giao phối]
    chr_best_fit_indx : [list]
        [Contains list of best fitness indices in chromosome population]
        [Chứa danh sách chỉ số thích nghi tốt nhất]
    pop : [numpy.ndarray]
        [numpy array of chromosome population to get best fitness chromosomes
         from last population]
        [Mảng chứa các pop có chỉ số thích nghi tốt nhất]
    Returns
    -------
    [numpy.ndarray]
        [numpy array of chromosome population undergone crossover]
    """

    crossover_pop = np.zeros((Config.pop_max, Config.chr_len))#mảng gồm các phần tử pop = 0 và có chr_len phần tử

    crossover_pop[0, :] = pop[chr_best_fit_indx[0], :]
    crossover_pop[1, :] = pop[chr_best_fit_indx[1], :]
    crossover_pop[2, :] = pop[chr_best_fit_indx[2], :]

    itr = 3

    while itr < Config.pop_max / 5:

        a = random.randint(0, Config.chr_len - 1)
        b = random.randint(0, Config.chr_len - 1)#gán giá trị cho a và b từ sinh số ngẫu nhiên

        partner_a = ranked_pop[a, :]
        partner_b = ranked_pop[b, :]#cắt 2 mảng
        joining_pt = random.randint(0, Config.chr_len - 1)#sinh số ngẫu nhiên từ 0 đến 15

        crossover_pop[itr, :joining_pt] = partner_a[:joining_pt]
        crossover_pop[itr+1, :joining_pt] = partner_b[:joining_pt]

        crossover_pop[itr, joining_pt:] = partner_b[joining_pt:]
        crossover_pop[itr+1, joining_pt:] = partner_a[joining_pt:]#nối mảng theo sinh số ngẫu nhiên

        itr += 2

    while itr < Config.pop_max:

        crossover_pop[itr] = ranked_pop[itr]
        itr += 1

    return crossover_pop