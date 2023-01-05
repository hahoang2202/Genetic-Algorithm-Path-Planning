from config import Config

import numpy as np


def ranking(chr_pop_fitness, pop):
    """
    Chức năng này gói gọn khả năng tạo thứ hạng của nhiễm sắc thể
    dân số dựa trên phương pháp lựa chọn bánh xe roulet.
    
    Parameters
    ----------
    chr_pop_fitness : [numpy.ndarray]
        [Chỉ số của nhiễm sắc thể]
    pop : [numpy.ndarray]
        [Quần thể nhiễm sắc thể sẽ trải qua quá trình xếp hạng]
    
    Returns
    -------
    [numpy.ndarray]
        [Population of ranked chromosomes]
    """

    chromo_prob = cal_prob(chr_pop_fitness=chr_pop_fitness)

    chromo_cum_prob = np.cumsum(chromo_prob, axis=0)

    chromo_rank = _ranking_based_on_roulet_wheel_selection(chr_cum_prob=chromo_cum_prob)

    chromo_pop_ranked = _generate_mating_pool(chr_rank=chromo_rank, pop=pop)

    return chromo_pop_ranked


def _generate_mating_pool(chr_rank, pop):
    """
    Chức năng này chịu trách nhiệm tạo nhóm giao phối sẽ trải qua
    lai ghép và đột biến ở giai đoạn tiếp theo.
    
    Parameters
    ----------
    chr_rank : [numpy.ndarray]
        [Bậc của các nhiễm sắc thể dựa trên phương pháp chọn lọc tối ưu]
    pop : [numpy.ndarray]
        [Quần thể nhiễm sắc thể từ đó nhóm giao phối sẽ được tạo ra]
    
    Returns
    -------
    [numpy.ndarray]
        [numpy arrayy of mating pool]
    """

    ranked_pop = np.zeros((1, np.shape(pop)[1]))

    for i in range(Config.pop_max):
        for j in range(int(chr_rank[i, 0])):
            if np.shape(ranked_pop)[0] == 1:
                ranked_pop = pop[i, :]
            else:
                ranked_pop = np.vstack((ranked_pop, pop[i, :]))
    return ranked_pop


def _ranking_based_on_roulet_wheel_selection(chr_cum_prob):
    """
    Chức năng này gói gọn khả năng thực hiện xếp hạng các nhiễm sắc thể
    dân số dựa trên phương pháp lựa chọn bánh xe
    
    Parameters
    ----------
    chr_cum_prob : [numpy.ndarray]
        [Xác suất tích lũy nhiễm sắc thể dựa trên thể lực của nhiễm sắc thể]
    
    Returns
    -------
    [numpy.ndarray]
        [bậc của các nhiễm sắc thể dựa trên phương pháp chọn lọc bánh xe roulets]
    """

    rand_array = np.random.rand(Config.pop_max) # mảng chứa các phần tử được sinh ngẫu nhiên
    no_of_times_chr_got_choosen = np.zeros((Config.pop_max, 1)) # mảng NST được chọn chứa các phần tử 0
    chr_rank = np.zeros((Config.pop_max, 1)) # mảng xếp hạng các NST được khởi tạo với các phần tử 0

    for i in range(Config.pop_max):
        k = 0
        while chr_cum_prob[k, 0] < rand_array[i]:
            k += 1
        no_of_times_chr_got_choosen[i, 0] = k

    for i in range(Config.pop_max):
        for j in range(Config.pop_max):
            if no_of_times_chr_got_choosen[j, 0] == i:
                chr_rank[i, 0] += 1

    return chr_rank


def cal_prob(chr_pop_fitness):
    """
    ----------
    chr_pop_fitness : [numpy.ndarray]
        [sự phù hợp của quần thể nhiễm sắc thể]
    
    Returns
    -------
    [numpy.ndarray]
        [Xác suất nhiễm sắc thể dựa trên chỉ số của nhiễm sắc thể]
    """

    chr_prob = np.zeros((Config.pop_max, 1))

    for i in range(Config.pop_max):

        chr_prob[i, 0] = chr_pop_fitness[i, 0] / \
            np.sum(chr_pop_fitness, keepdims=True)

    return chr_prob
