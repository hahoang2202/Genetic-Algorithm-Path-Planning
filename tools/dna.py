from config import Config

import numpy as np
import random


def dna(chr_pop_fitness, ranked_population, chr_best_fitness_index, last_pop):
    """
    Chức năng này gói gọn chức năng liên quan đến dna như chéo
    và đột biến.
    ----------
    chr_pop_fitness : [numpy.ndarray]
        [Chứa các giá trị thể chất của quần thể nhiễm sắc thể]

    ranked_population : [numpy.ndarray]
        [Chứa mảng numpy của quần thể nhiễm sắc thể được xếp hạng]

    chr_best_fitness_index : [list]
        [Chứa danh sách các chỉ số tốt nhất trong quần thể nhiễm sắc thể]

    last_pop : [numpy.ndarray]
        [Chứa mảng của dân số cuối cùng]
    Returns
    -------
    [numpy.ndarray]
        [numpy array of chromosome with have gone through random crossover and mutation]
    """

    chromo_crossover_pop = giao_phối(
        ranked_pop=ranked_population, DS_chỉ_số_cao_nhất=chr_best_fitness_index,
        pop=last_pop) #quá trình giao phối ADN với các ADN có chỉ số cao nhất

    chromo_crossover_mutated_pop = đột_biến(pop=chromo_crossover_pop)#giao phối sau khi đột biến

    return chromo_crossover_mutated_pop


def đột_biến(pop):#tạo đột biến
    """
    Chức năng này chịu trách nhiệm xử lý đột biến trong quần thể nhiễm sắc thể.
    ----------
    pop : [numpy.ndarray]
        [Mảng chứa các nhiễm sắc thể sẽ trải qua đột biến]
    Returns
    -------
    [numpy.ndarray]
        [Mảng chứa NST bị đột biến]
    """

    NST_đột_biến = np.array(pop, copy=True)#mảng NST đột biến

    itr = 3
    while itr < Config.pop_max: 
        for k in range(Config.chr_len):
            c = random.random()# sinh số ngẫu nhiên
            if c < Config.mutation_rate and k != 0:
                NST_đột_biến[itr, k] = random.randint(1, Config.npts - 2)
                '''
                Nếu tỉ lệ đột biến bé hơn 0.01 thì thay thế NST cũ bằng
                NST mới trừ 2 điểm đầu và cuối
                '''
            else:
                pass
        itr += 1
    return NST_đột_biến


def giao_phối(ranked_pop, DS_chỉ_số_cao_nhất, pop):#giao phối
    """
    Chức năng này chịu trách nhiệm xử lý trao đổi chéo trong quần thể nhiễm sắc thể
    ----------
    ranked_pop : [numpy.ndarray]
        [Mảng chứa các pop sẽ giao phối]
    DS_chỉ_số_cao_nhất : [list]
        [Chứa danh sách chỉ số thích nghi tốt nhất]
    pop : [numpy.ndarray]
        [Mảng chứa các pop có chỉ số thích nghi tốt nhất]
    Returns
    -------
    [numpy.ndarray]
        [numpy array of chromosome population undergone crossover]
    """

    NST_giao_phốt = np.zeros((Config.pop_max, Config.chr_len))#mảng gồm các phần tử pop = 0 và có chr_len phần tử

    NST_giao_phốt[0, :] = pop[DS_chỉ_số_cao_nhất[0], :]
    NST_giao_phốt[1, :] = pop[DS_chỉ_số_cao_nhất[1], :]
    NST_giao_phốt[2, :] = pop[DS_chỉ_số_cao_nhất[2], :]

    itr = 3

    while itr < Config.pop_max / 5:

        a = random.randint(0, Config.chr_len - 1)
        b = random.randint(0, Config.chr_len - 1)#gán giá trị cho a và b từ sinh số ngẫu nhiên

        partner_a = ranked_pop[a, :]
        partner_b = ranked_pop[b, :]#cắt 2 mảng
        joining_pt = random.randint(0, Config.chr_len - 1)#sinh số ngẫu nhiên từ 0 đến 15

        NST_giao_phốt[itr, :joining_pt] = partner_a[:joining_pt]
        NST_giao_phốt[itr+1, :joining_pt] = partner_b[:joining_pt]

        NST_giao_phốt[itr, joining_pt:] = partner_b[joining_pt:]
        NST_giao_phốt[itr+1, joining_pt:] = partner_a[joining_pt:]#nối mảng theo sinh số ngẫu nhiên

        itr += 2

    while itr < Config.pop_max:

        NST_giao_phốt[itr] = ranked_pop[itr]
        itr += 1

    return NST_giao_phốt