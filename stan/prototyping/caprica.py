from collections import namedtuple
import ipywidgets as widgets
import pystan
import glue
import json
import pandas as pd
import numpy as np
import pickle
import bqplot as bq

try:
    with open("predict.pkl", "rb") as f:
        sm_predict = pickle.load(f)
except FileNotFoundError:
    sm_predict = pystan.StanModel("predict.stan")
    with open("predict.pkl", "wb") as f:
        pickle.dump(sm_predict, f)
try:
    with open("infer.pkl", "rb") as f:
        sm_infer = pickle.load(f)
except FileNotFoundError:
    sm_infer = pystan.StanModel("infer.stan")
    with open("infer.pkl", "wb") as f:
        pickle.dump(sm_infer, f)

class LegendWidget:
    """A legend Widget using a horizontal bar chart
    
    Adapted from DougRzz on GitHub:
    https://github.com/bloomberg/bqplot/issues/601

    marks: line marks from a bqplot figure.  
    
    These line marks must have legend labels 
    (in line mark, remove other legend by using this: display_legend = False)
    e.g. >>> legend = legendWidget(fig.marks) 
    
    BQplot module imported as bq (import bqplot as bq)

    """
    def __init__(self, marks):
        """Return a new Legend object."""
        y_ord = bq.OrdinalScale()
        x_sc = bq.LinearScale()
        
        legendLabels = []
        colours = []
        markLineNums = [] # record number of lines per mark
        # marks.reverse()
        for mark in marks:            
            legendLabels += mark.labels
            colours += mark.colors[:len(mark.labels)]
            markLineNums.append(len(mark.labels))  

        bar = bq.Bars(
            y=[1]*len(legendLabels) , # all bars have a amplitude of 1
            x=legendLabels, 
            scales={'y': x_sc, 'x': y_ord},
            colors=colours ,
            padding = 0.6,
            orientation='horizontal',
            stroke = 'white'  #remove the black border around the bar
            )
        
        ax_y = bq.Axis(scale=y_ord, orientation="vertical")
        ax_x = bq.Axis(scale=x_sc)
        ax_x.visible = False
        margin = dict(top=40, bottom=0, left=110, right=5)
        barFig = bq.Figure(marks=[bar], axes=[ax_y, ax_x], fig_margin=margin)
        
        # Variable height depending on number of bars in legend
        barFig.layout.height = str(45 + 20 * len(legendLabels)) + 'px'
        barFig.layout.width = '170px'

        barFig.min_aspect_ratio = 0.000000000001 # effectively remove aspect ratio constraint
        barFig.max_aspect_ratio = 999999999999999 # effectively remove aspect ratio constraint
        barFig.background_style = {'fill': 'White'}   
                    
        self.fig = barFig
        self.bar = bar
        self.colours = colours
        self.markLineNums = markLineNums
   
PredictedOutcome = namedtuple("PredictedOutcome", ["t", "low", "med", "high"])

class PredictiveInteractive(widgets.VBox):
    
    def __init__(self, bm_weights=None, infer_csv="inferred.csv", infer_json="infer_input.json",
                 time_samples=np.arange(0, 7, 1)):
        self._infer_csv = infer_csv
        self._infer_json = infer_json
        self._time_samples = time_samples
        with open(infer_json) as f:
            data = json.load(f)
        self.N_bm = data["N_bm"]
        self.N_tx = data["N_tx"]
        self.outcomes = None
        if bm_weights is None:
            bm_weights = np.ones(self.N_bm) / self.N_bm
        self._bm_weights = bm_weights
        
        # widgets!
        biomarker_boxes = widgets.VBox([widgets.Checkbox(description=f"Biomarker {i + 1}", value=False)
                                        for i in range(self.N_bm)])
        self.biomarker_boxes = biomarker_boxes
        
        randomize_button = widgets.Button(description="Random Biomarkers")
        randomize_button.on_click(self.randomize_biomarkers)
        plot_button = widgets.Button(description="Plot")
        plot_button.style.button_color = "lightgreen"
        plot_button.on_click(self.plot_outcomes)
        buttons = widgets.VBox([randomize_button, plot_button])
        input_ui = widgets.HBox([buttons, biomarker_boxes])
        
        xscale = bq.LinearScale()
        yscale = bq.LinearScale()
        scales = dict(x=xscale, y=yscale)
        marks = []
        t = time_samples
        self._colors = ["blue", "red", "green", "orange", "purple"]
        for tx in range(self.N_tx):
            c = self._colors[tx]
            line = bq.Lines(x=t, y=np.zeros_like(t), scales=scales,
                            labels=[f"Treatment {tx + 1}"],
                            colors=[c],
                            opacities=[],
                            display_legend=False)
            marks.append(line)
            line = bq.Lines(x=t,
                            y=(-np.ones_like(t), np.ones_like(t)),
                            scales=scales,
                            fill="between",
                            fill_opacities=[0.2],
                            opacities=[0.0, 0.0],
                            colors=[c],
                            display_legend=False)
            marks.append(line)

        fig = bq.Figure(marks=marks,
                        title="Predicted Outcomes",
                        legend_location="bottom-left",
                        axes=[bq.Axis(scale=xscale, label="Time (months)"),
                              bq.Axis(scale=yscale, label="Change in Response",
                                      orientation="vertical")])

        legend = LegendWidget(fig.marks)
        legend.fig.marks[0].on_element_click(self.change_opacity)
        self._legend = legend
        self._fig = fig
        figures = widgets.HBox([fig, legend.fig])
        super().__init__([input_ui, figures])
        
    def randomize_biomarkers(self, useless_input):
        for bm in range(self.N_bm):
            p = self._bm_weights[bm]
            value = bool(np.random.choice([False, True], p=[1 - p, p]))
            self.biomarker_boxes.children[bm].value = value
            
    def change_opacity(self, bar, target):
        """Enable legend interactivity. 
        Use in conjunction with class legendWidget(object) 
        Click on legend bar to toggle opacity of all other lines

        """
        lineFig = self._fig 
        legendFig = self._legend.fig

        opacity = 0.1   # set opacity of non selected lines here
        bar_index = target['data']['index']
        # index = self.N_tx - bar_index # was reversed
        index = bar_index
        
        if bar.opacities == [] or bar.opacities[index] == opacity:
            # highlight indexed line
            bar.opacities=[opacity]*index + [1] + [opacity]*(len(bar.x) - index - 1)        
            for i in range(self.N_tx):
                med_mark = lineFig.marks[2 * i + 0]
                fill_mark = lineFig.marks[2 * i + 1]
                if i == index:
                    med_mark.opacities = []
                    fill_mark.fill_opacities = [0.2]
                    fill_mark.opacities = [0.5, 0.5]
                else:
                    med_mark.opacities = [0.1]
                    fill_mark.fill_opacities = [0.05]
                    fill_mark.opacities = [0.0, 0.0]
        else:
            # reset lines
            bar.opacities = []
            for i in range(self.N_tx):
                med_mark = lineFig.marks[2 * i + 0]
                fill_mark = lineFig.marks[2 * i + 1]
                med_mark.opacities = []
                fill_mark.fill_opacities = [0.2]
                fill_mark.opacities = [0.0, 0.0]

                
    def generate_predictions(self):
        bm_indicators = [np.array([bm_box.value for bm_box in self.biomarker_boxes.children]).astype(int).tolist()]
        data = glue.predict_from_infer(self._infer_csv,
                                       self._infer_json,
                                       outfile=None,
                                       bm_indicators=bm_indicators,
                                       time_samples=self._time_samples)
        fit = sm_predict.sampling(data=data, chains=1, iter=1000, warmup=500)
        pt_index = np.array(data["pt_index"])
        N_res = np.array(data["N_res"])
        N_bm = np.array(data["N_bm"])
        N_tx = np.array(data["N_tx"])
        bm_indicators = np.array(data["bm_indicators"])
        tx_indicators = np.array(data["tx_indicators"])
        t_res = np.array(data["t_res"])
        self.outcomes = {}
        for tx in range(1, 1 + N_tx):
            idx = np.arange(N_res)[pt_index == tx]
            t = t_res[idx]
            pars = [f"y_res[{i + 1}]" for i in idx]
            samples = fit.extract(pars=pars)
            y_res = np.array(list(samples.values()))
            low, med, high = np.percentile(y_res, axis=1, q=[16, 50, 84])
            self.outcomes[tx] = PredictedOutcome(t, low, med, high)
            
    def plot_outcomes(self, *args, **kwargs):
        self.generate_predictions()
        for i, (tx, outcome) in enumerate(self.outcomes.items()):
            line = self._fig.marks[2 * i + 0]
            line.x = outcome.t            
            line.y = outcome.med
            line = self._fig.marks[2 * i + 1]
            line.x = outcome.t
            line.y = (outcome.low, outcome.high)
        return self._fig

    
