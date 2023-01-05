from config import Config

import numpy as np
import math as ma
import random


def population():
    """
    Chức năng này gói gọn khả năng khởi tạo quần thể nhiễm sắc thể.
    
    Returns
    -------
    [numpy.ndarray]
        [Population of chromosomes]
    """
    
    link = Config.define_links()
    link_fit = _link_distance(link)
    link_prob = _link_prob(link_fit)
    link_cum_prob = np.cumsum(link_prob, axis=1)
    initial_pop = _create_pop(link_cum_prob=link_cum_prob)

    return initial_pop


def _link_distance(link):
    """
    Chức năng này chịu trách nhiệm tính toán khoảng cách giữa 2 điểm có thể liên kết
    
    Parameters
    ----------
    link : [numpy.ndarray]
        [Liên kết giữa các điểm]
    
    Returns
    -------
    [numpy.ndarray]
        [numpy array of distance b/w links]
    """

    link_dist = np.zeros((np.shape(link)[0], np.shape(link)[1]-1))

    for i in range(np.shape(link)[0]):

        for j in range((np.shape(link)[1])-1):

            if link[i][j] > -0.1 and link[i][j+1] > -0.1:

                link_dist[i][j] = calculate_distance (
                    pt_1=Config.path_points[int(link[i][j])],
                    pt_2=Config.path_points[int(link[i][j+1])]
                )
                    
    return link_dist


def _link_prob(link_fit):
    """
    Chức năng này tính toán xác suất của các liên kết.
    
    Parameters
    ----------
    link_fit : [numpy.ndarray]
        [mảng numpy của các kết nối liên kết tập chỉ số dựa trên khoảng cách]
    
    Returns
    -------
    [numpy.ndarray]
        [mảng numpy của xác suất liên kết dựa trên độ phù hợp của liên kết]
    """

    link_prob = np.zeros((np.shape(link_fit)[0], np.shape(link_fit)[1]))

    for i in range(np.shape(link_fit)[0]):

        for j in range(np.shape(link_fit)[1]):

            link_prob[i][j] = link_fit[i][j]/np.sum(link_fit[i], keepdims=True)

    return link_prob


def _create_pop(link_cum_prob):
    """
    Chức năng này chịu trách nhiệm tạo ra tập hợp nhiễm sắc thể dựa trên sự kết nối
    liên kết b/w.
    
    Parameters
    ----------
    link_cum_prob : [numpy.ndarray]
        [mảng numpy của xác suất tích lũy liên kết dựa trên độ phù hợp của liên kết]
    
    Returns
    -------
    [numpy.ndarray]
        [mảng numpy của tập thể nhiễm sắc thể dựa trên kết nối liên kết b/w]
    """
    pop = np.zeros((Config.pop_max, Config.chr_len))
    pop[:, 0] = Config.start_index
    pop[:, Config.chr_len - 1] = Config.end_index

    link = Config.define_links()

    for k in range(Config.pop_max):
        i = Config.start_index
        j = Config.start_index + 1
        while j < Config.chr_len:
            i = int(i)
            if j > 0 and j < (Config.chr_len - 1):
                random_val = random.random()
                if random_val < link_cum_prob[i][0]:
                    pop[k][j] = link[i][1]
                    i = link[i][1]
                    if _both_equ(i, Config.end_index):
                        while j < (Config.chr_len - 1):
                            pop[k][j+1] = Config.end_index
                            j += 1
                elif random_val < link_cum_prob[i][1]:
                    pop[k][j] = link[i][2]
                    i = link[i][2]
                    if _both_equ(i, Config.end_index):
                        while j < (Config.chr_len - 1):
                            pop[k][j+1] = Config.end_index
                            j += 1
                elif random_val < link_cum_prob[i][2]:
                    pop[k][j] = link[i][3]
                    i = link[i][3]
                    if _both_equ(i, Config.end_index):
                        while j < (Config.chr_len - 1):
                            pop[k][j+1] = Config.end_index
                            j += 1
                elif random_val < link_cum_prob[i][3]:
                    pop[k][j] = link[i][4]
                    i = link[i][4]
                    if _both_equ(i, Config.end_index):
                        while j < (Config.chr_len - 1):
                            pop[k][j+1] = Config.end_index
                            j += 1
            j += 1
    return pop


def _both_equ(element_1, element_2):
    """
    Hàm này chịu trách nhiệm tìm xem cả hai phần tử có bằng nhau hay không.
    
    Parameters
    ----------
    element_1 : [Int]
        [Yếu tố đầu tiên để so sánh]
    element_2 : [Int]
        [Yếu tố thứ hai để so sánh]
    
    Returns
    -------
    [Bool]
        [Đúng hay Sai? Cả hai yếu tố có bằng nhau hay không?]
    """

    return True if int(element_1) == int(element_2) else False


def calculate_distance(pt_1, pt_2):
    """
    Chức năng này gói gọn khả năng tính toán khoảng cách b/w hai điểm.
    
    Parameters
    ----------
    pt_1 : [Float]
        [điểm 1 để tính khoảng cách]
    pt_2 : [Float]
        [điểm 2 để tính khoảng cách]
    
    Returns
    -------
    [float]
        [Khoảng cách 2 điểm]
    """
    return ma.sqrt(ma.pow((pt_1[0]-pt_2[0]), 2) + ma.pow((pt_1[1]-pt_2[1]), 2))# độ dài vecto

