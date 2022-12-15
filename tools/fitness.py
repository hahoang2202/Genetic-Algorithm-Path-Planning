"""
Script provides functionality related to fitness of chromosomes based
on total distance and connection b.w two consecutive nodes of a chromosome

Author: Yasim Ahmad(yaaximus)

Email: yasim.ahmed63@yahoo.com
"""

from config import Config
from tools.population import calculate_distance

import numpy as np


def fitness(chr_pop):
    """
    This function is responsible for calculating fitness of chromosomes.
    Chức năng này chịu trách nhiệm tính toán mức độ phù hợp của NST.
    
    Parameters
    ----------
    chr_pop : [numpy.ndarray]
        [Population of chromosomes whose fitness is to be calculated]
        Các NST có thể tính toán
    
    Returns
    -------
    [numpy.ndarray]
        [(1)Population of chromosomes whose fitness is calculated, 
         (2)List of indices of best fitness chromosomes]
    """

    chromo_pts_consec_dist = chr_pts_consecutive_dist(pop=chr_pop)

    chromo_fit_based_dist = chr_fit_based_dist(chr_pts_consec_dist=chromo_pts_consec_dist)

    chromo_conn = chr_conn(chr_pop=chr_pop)

    chromo_fit_based_conn = chr_fit_based_conn(chr_conn=chromo_conn)

    chromo_fit = chr_fit(chr_fit_based_dist=chromo_fit_based_dist,
        chr_fit_based_conn=chromo_fit_based_conn)

    chromo_best_fit_index = chr_best_fit_ind(chr_fit=chromo_fit)

    return chromo_fit, chromo_best_fit_index


def chr_best_fit_ind(chr_fit):
    """
    This function is responsible for finding best fitness chromosome indices.
    Chức năng này chịu trách nhiệm tìm kiếm các chỉ số thích nghi tốt nhất NST.
    
    Parameters
    ----------
    chr_fit : [numpy.ndarray]
        [numpy array of chromosome population fitness]
        Mảng chỉ số thích nghi
    
    Returns
    -------
    [list]
        [list of best fitness chromosome indices]
        Các NST có chỉ số tốt nhất
    """

    temp_chr_fit = np.array(chr_fit, copy=True)# mảng các chỉ số

    chr_best_fit_index = []

    while len(chr_best_fit_index) < 3:

        y = np.where(temp_chr_fit == np.amax(temp_chr_fit))[0]

        for i in range(len(y)):
            chr_best_fit_index.append(int(y[i]))# thêm vào cuối mảng

        for i in chr_best_fit_index:
            temp_chr_fit[i][0] = 0

    return chr_best_fit_index


def chr_fit(chr_fit_based_dist, chr_fit_based_conn):
    """
    This function is responsible for calculating fitness of chromosome population
    based on total distance of individual chromosome, and links b/w path points
    of individual chromosome.
    
    Tính toán mức độ phù hợp của NST giữa trên khoảng cách và khả năng kết nối giữa 2 NST
    
    Parameters
    ----------
    chr_fit_based_dist : [numpy.ndarray]
        [numpy array of chromosome fitness based on total distance]
        Mảng chỉ số khoảng cách
    chr_fit_based_conn : [numpy.ndarray]
        [numpy array of chromosome fitness based on links b/w path 
        points of individual chromosome]
        Mảng khả năng kết nối giữa 2 NST (2 điểm)
    
    Returns
    -------
    [numpy.ndarray]
        [final fitness of chromosome population]
    """

    chr_fit = np.zeros((Config.pop_max, 1))

    for i in range(Config.pop_max):

        chr_fit[i][0] = chr_fit_based_dist[i][0] + chr_fit_based_conn[i][0]

    return chr_fit


def chr_fit_based_conn(chr_conn):
    """
    This function is responsible for calculating fitness of chromosome population
    based on number of connections b/w path points of an individual chromosome
    Tính toán mức độ phù hợp dựa trên số lượng kết nối từ 1 NST đến các NST khác
    
    Parameters
    ----------
    chr_conn : [numpy.ndarray]
        [numpy array of number of connection b/w path points of an individual chromosome]
        Mảng chứa số lượng điểm có thể kết nối từ 1 điểm
    
    Returns
    -------
    [numpy.ndarray]
        [numpy array of chromosome population fitness based on connections]
    """

    chr_conn_fit = np.zeros((Config.pop_max, 1))

    for i in range(Config.pop_max):
    #chỉ số thích nghi được tính dựa trên tỉ số SL kết nối / SL kết nối tối đa
        chr_conn_fit[i][0] = chr_conn[i][0] / ( Config.chr_len - 1 )

    return chr_conn_fit


def chr_conn(chr_pop):
    """
    This function is responsible for finding number of connections b/w path points of
    a individual chromosome.
    Hàm tìm kiếm số đường dẫn từ 1 điểm
    
    Parameters
    ----------
    chr_pop : [numpy.ndarray]
        [Population of chromosomes whose number of connections are to be calculated]
        Các NST cần tính toán số lượng đường dẫn
    
    Returns 
    -------
    [numpy.ndarray]
        [numpy array of number of connection b/w path points of an individual chromosome]
        Số lượng đường dẫn của 1 điểm
    """

    link = Config.define_links()
    chr_conn = np.zeros((Config.pop_max, 1))

    for i in range(Config.pop_max):
        for j in range(Config.chr_len-1):
            a = int(chr_pop[i][j])
            b = int(chr_pop[i][j+1])
            for k in range(np.shape(link)[1]):# trả về kích thước của mảng
                if link[a, k] == b:# thỏa mãn có thể kết nối
                    chr_conn[i][0] += 1

    return chr_conn


def chr_fit_based_dist(chr_pts_consec_dist):
    """
    This function is responsible for calculating chromosome fitness based on total
    distance of individual chromosome.
    Tính toán dựa trên tổng khoảng cách
    
    Parameters
    ----------
    chr_pts_consec_dist : [numpy.ndarray]
        [numpy array of individual chromosome total distance]
        Chứa tổng khoảng cách của 1 điểm
    
    Returns
    -------
    [numpy.ndarray]
        [numpy array of individual chromosome fitness based on total distance]
    """

    chr_pop_fit_based_dist = np.zeros((Config.pop_max, 1))

    for i in range(Config.pop_max):

        chr_pop_fit_based_dist[i][0] = 10.0 * \
            (1.0 / np.sum(chr_pts_consec_dist[i], keepdims=True))

    return chr_pop_fit_based_dist


def chr_pts_consecutive_dist(pop):
    """
    This function is responsible for calculating total distance of individual
    chromosome in population of chromosomes.
    Tính tổng khoảng cách của NST
    
    Parameters
    ----------
    pop : [numpy.ndarray]
        [Population of chromosomes whose total distance is to be calculated]
        Các NST
    
    Returns
    -------
    [numpy.ndarray]
        [numpy array of individual chromosome total distance]
    """

    chr_pop_dist = np.zeros((Config.pop_max, Config.chr_len-1))

    for i in range(Config.pop_max):

        for j in range(Config.chr_len-1):

            chr_pop_dist[i][j] = calculate_distance(
                pt_1=Config.path_points[int(pop[i][j+1])],
                pt_2=Config.path_points[int(pop[i][j])]) #Hàm tính khoảng cách giữa 2 điểm

    return chr_pop_dist
