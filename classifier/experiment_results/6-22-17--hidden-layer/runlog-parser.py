import ast, csv
# import matplotlib.pyplot as plt

class DataLine:
    def __init__(self, d):
        self.data = ast.literal_eval(d)

    def write_line(self, writer, fields):
        vals = [self.data[field] for field in fields]
        writer.writerow(vals)

class DataCollection:
    def __init__(self, path_to_runlog):
        self.datafile = open(path_to_runlog, 'r')
        self.lines = self.datafile.readlines()
        self.data = [DataLine(line) for line in self.lines]

class CSVWriter:
    def __init__(self, path_to_runlog, path_to_output, fields):
        self.data_collection = DataCollection(path_to_runlog)
        self.csvfile = open(path_to_output, 'w')
        self.writer = csv.writer(self.csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
        self.fields = fields

    def build_csv(self):
        self.write_header()
        self.write_data()
        self.csvfile.close()

    def write_header(self):
        self.writer.writerow(self.fields)

    def write_data(self):
        for dl in self.data_collection.data:
            dl.write_line(self.writer, self.fields)

#under development by nathaniel...
# class Plotter:
#     def __init__(self, path_to_runlog):
#         self.data_collection = DataCollection(path_to_runlog)


#ALLOWED FIELDS
#'normal_training_divisor', 'training_mode', 'run_label', 'GD_step_size', 'initial_weights_and_bias_mode',
#'random_background_shift', 'ts',  'results', 'random_training_shift', 'trial_n', 
#'stretch_training_divisor', 'dim_of_input', 'random_background_multiplier', 'stretch_training_shift', 'training_constant', 
#'background_constant', 'random_training_multiplier', 'mean_results', 'test_fraction', 'background_mode'
#'testing_length', 'delta_k_training_inclusion_abs_threshold', 'delta_k_testing_inclusion_abs_threshold', 'training_length'
fields = ['training_mode', 'mean_results', 'results'] 
csvw = CSVWriter('runlog.dbl', 'data.csv', fields)
csvw.build_csv()
