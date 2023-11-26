#  
#   main program :
#   life origin simulator
#
from lifeorig.parser import parser
from lifeorig.read_input import p
from lifeorig.chromosomes_popul import chromosomes_pop_class
from lifeorig.selection_methods import tournament_selection

# log.info(print_header())

args = parser.parse_args()
p.read_input_json(args.json_input)

# set up GA

chromosomes = chromosomes_pop_class()
chromosomes.set_initial_population()
chromosomes.compute_scores()
chromosomes.print_avg_score()

for iter in range(p.n_max_iter):
    selected = [tournament_selection(chromosomes) for _ in range(chromosomes.n_pop)]
    chromosomes.produce_next_gen(selected)
    chromosomes.compute_scores()
    chromosomes.print_avg_score()