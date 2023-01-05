from tools.population import population
from tools.fitness import fitness
from tools.ranking import ranking
from tools.dna import dna
from tools.draw_plot import show_plot

from config import Config

def main():
    """
    Chức năng này gói gọn khả năng khởi tạo quần thể nhiễm sắc thể
    và sau đó tiếp tục tính toán chỉ số thích nghi, tạo xếp hạng, thực hiện chéo và
    đột biến & lặp lại các bước trên cho đến khi không đáp ứng tiêu chí dừng đã xác định
    """

    chr_population = population()#nhóm các NST

    chr_pop_fitness, chr_best_fitness_index = fitness(chr_pop=chr_population) # tính toán chỉ sô thích nghi

    chr_ranked_population = ranking(chr_pop_fitness=chr_pop_fitness, pop=chr_population) # xếp hạng các chỉ số theo nhóm

    chr_crossover_mutated_population = dna(chr_pop_fitness=chr_pop_fitness,
        ranked_population=chr_ranked_population, chr_best_fitness_index=
        chr_best_fitness_index, last_pop=chr_population) # tiến hành đột biến các NST từ dưới lên

    show_plot(best_chromosome=chr_crossover_mutated_population[0]) # tiến hành vẽ đồ thị chứa đường dẫn tối ưu

    while not Config.stop_generation:# tiêu chí dừng

        prev_best_fit = chr_pop_fitness[chr_best_fitness_index[0], 0]

        chr_pop_fitness, chr_best_fitness_index = fitness(
            chr_pop=chr_crossover_mutated_population)

        chr_ranked_population = ranking(chr_pop_fitness=chr_pop_fitness, 
            pop=chr_crossover_mutated_population)

        chr_crossover_mutated_population = dna(chr_pop_fitness=chr_pop_fitness,
            ranked_population=chr_ranked_population, chr_best_fitness_index=
            chr_best_fitness_index, last_pop=chr_crossover_mutated_population)

        if prev_best_fit == chr_pop_fitness[chr_best_fitness_index[0], 0]:
            Config.stop_criteria += 1 # tăng tiêu chuẩn dừng lên 1
        else:
            Config.stop_criteria = 0 # đặt lại về 0

        if Config.stop_criteria >= 5: # khi tiêu chí đạt lớn hơn 4 thì dừng
            Config.stop_generation = True

        print("Best chromosome is:", chr_crossover_mutated_population[chr_best_fitness_index[0]])

        show_plot(best_chromosome=chr_crossover_mutated_population[0])
        Config.generations += 1

    show_plot(best_chromosome=chr_crossover_mutated_population[0], inf_time=True)

if __name__ == '__main__':

    main()
