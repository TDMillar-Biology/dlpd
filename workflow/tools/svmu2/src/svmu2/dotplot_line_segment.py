''' TDMILLAR MAY 2024 A dot plot line segment class for storing line segments that bridge gaps between 'real' lines on a dot plot
used to determine what the gap represents (ie insertion / deletion)
'''
from matplotlib import pyplot as plt
from math import sqrt
class DotPlotLineSegment:
    def __init__(self, chrom, reference_start, reference_end, query_start, query_end, sv_type, theta, index):
        #CHROM POS      ID         REF   ALT    QUAL  FILTER   INFO
        self.chrom = chrom
        self.reference_start = int(reference_start)
        self.reference_end = int(reference_end)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.sv_type = sv_type
        self.x_mid = (self.reference_start + self.reference_end) / 2
        self.y_mid = (self.query_start + self.query_end) / 2
        self.shares_domain = []
        self.theta = theta
        self.length = sqrt((self.reference_end - self.reference_start)**2 + (self.query_end - self.query_start)**2)
        self.range_partners = []
        self.domain_partners = []
        self.index = index
        self.event_ID = f'{chrom}_{index:04}'
        try:
            self.slope = (self.query_end - self.query_start) / (self.reference_end - self.reference_start)
        except ZeroDivisionError:
            self.slope = float('inf')
        
    def __str__(self):
        #return f'{self.chrom}\t{self.reference_start}\tN\t{self.sv_type}\taddqual\taddfilter\taddinfo'
        return f'SV Type: {self.sv_type}\n Theta: {self.theta}\n ref: {(self.reference_start, self.reference_end)} \n query: {(self.query_start, self.query_end)}\n'
    
    def plot(self, color = 'orange', text = False, components = False):
        plt.plot([self.reference_start, self.reference_end], [self.query_start, self.query_end], color=color)
        if components:
            plt.plot([self.reference_start, self.reference_start], [self.query_start, self.query_end], color = 'blue') ## plot the x component
            plt.plot([self.reference_start, self.reference_end], [self.query_start, self.query_start], color = 'blue') ## plot the y component
        #plt.text(self.x_mid, self.y_mid, f'SV Type: {self.sv_type}')
        if text:
            plt.text(self.x_mid, self.y_mid, f'SV Type: {self.sv_type}\n Theta: {self.theta}', fontsize=5, color="black")
            plt.text(self.x_mid, self.y_mid, f'SV Type: {self.sv_type}\n Theta: {self.theta}\n ref: {(self.reference_start, self.reference_end)} \n query: {(self.query_start, self.query_end)}', fontsize=8, color="black")

