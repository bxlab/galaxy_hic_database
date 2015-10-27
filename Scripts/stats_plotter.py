#!/usr/bin/env python

import sys
import os

import h5py
import numpy


def main():
    data_fname, out_fname = sys.argv[1:3]
    data = h5py.File(data_fname, 'r')
    template = """
    <!DOCTYPE html>
    <html lang="en-US">

    <link href="http://cdn.pydata.org/bokeh/release/bokeh-0.9.0.min.css"
        rel="stylesheet" type="text/css">
    <script src="http://cdn.pydata.org/bokeh/release/bokeh-0.9.0.min.js"></script>
    <script type="text/javascript" src="https://www.google.com/jsapi"></script>
    <head>
    <title>
    """ + (data_fname.split('/')[-1].split('.')[0]) + """
     </title>
    <style>
    h4, h2 {
        padding-top: 10px;
    }
    th {
    text-align: center;
        border-bottom: 1px solid black;
    }
    th:nth-child(1) {
        text-align: right;
    }
    td{
        border-bottom: 1px solid black;
        border-collapse: collapse;
        text-align: center;
        padding-left: 8px;
        padding-right: 8px;
        padding-bottom: 5px;
        padding-top: 5px;
    }
    td:nth-child(1) {
        text-align: right;
        font-weight: bold;
    }
    tr:nth-child(1) {
        font-weight: bold;
        text-align: right;
    }
    </style>
    </head>
    <body>

    """
    for arg in sys.argv[3:]:
        name, value = arg.split('=')
        template += "<p><b>%s:</b> %s</p>\n" % (name.replace('_', ' '), value.strip('"').replace('_', ' '))
    template += "<p><b>Sample ID:</b> %s</p>\n" % (data_fname.split('/')[-1].split('.')[0])
    template += load_mapping_stats(data)
    template += load_hcd_stats(data)
    best_insert, html = load_insert_distribution(data)
    template += html
    template += load_cis_interaction_distribution(data)
    template += load_trans_interaction_distribution(data)
    template += """<h3>Best insert size</h3>\n%i\n""" % int(round(best_insert))    
    template += """</body>\n</html>"""
    output = open(out_fname, 'w')
    output.write(template)
    output.close()
    print "%s\t%i" % (data_fname.split('/')[-1].split('.')[0], int(round(best_insert)))

def load_mapping_stats(data):
    valid = False
    stats = []
    if 'raw_filelist' in data['/'].attrs.keys():
        files = data['/'].attrs['raw_filelist'].split(',')
        if files[0].split('.')[-1] == 'raw':
            valid = True
        if valid:
            for fname in files:
                if not os.path.exists(fname.replace('raw', 'stats')):
                    valid = False
    if not valid:
        return ""
    x_labels = ['File']
    y_labels = []
    stats = []
    for i, fname in enumerate(files):
        stats.append([])
        x_labels.append(fname.split('/')[-1].split('.')[0])
        for line in open(fname.replace('raw', 'stats')):
            temp = line.rstrip('\n').split(': ')
            if temp[0].count(' quality ') > 0:
                continue
            if i == 0:
                y_labels.append(temp[0])
            temp = temp[1].split(" (")
            if len(temp) > 1:
                stats[-1].append("%s (%s" % (int2str(int(temp[0])), temp[1]))
            else:
                stats[-1].append(int2str(int(temp[0])))
    html = "<h4>Mapping Statistics</h4>\n"
    i = 1
    while i < len(x_labels):
        html += make_table([x_labels[0]] + x_labels[i:min(len(x_labels), i + 5)], y_labels,
                           stats[(i - 1):(min(len(stats), i + 4))]) + "<p></p>\n"
        i += 5
    return html + "<p></p>"

def load_hcd_stats(data):
    stats_table = data['stats'][...]
    label_dict = {
        'total_reads': ['Total reads', 0],
        'chr_not_in_fends': ['Chromosome not in fends', 1],
        'pcr_duplicates': ['PCR duplicates', 2],
        'out_of_bounds': ['Outide RE cutsite range', 3],
        'insert_size': ['Insert larger than cutoff', 4],
        'same_fragment': ['Fragment circularization', 5],
        'failed_cut': ['Possible failed RE digest', 6],
        'valid_cis_reads': ['Valid cis reads', 7],
        'valid_trans_reads': ['Valid trans reads', 8],
        'valid_cis_pairs': ['Unique cis fend pairs', 9],
        'valid_trans_pairs': ['Unique trans fend pairs', 10],
    }
    y_labels = []
    stats = [[]]
    for i in range(stats_table.shape[0]):
        y_labels.append('')
        stats[0].append(0)
    for name, count in stats_table:
        label, pos = label_dict[name]
        y_labels[pos] = label
        stats[0][pos] = count
    for i in range(1, 9):
        stats[0][i] = "%s (%0.2f%%)" % (int2str(stats[0][i]), 100.0 * stats[0][i] / stats[0][0])
    for i in [0, 9, 10]:
        stats[0][i] = int2str(stats[0][i])
    html = "<h4>HiC Read Statistics</h4>\n" + make_table([], y_labels, stats)
    return html

def make_table(x_labels, y_labels, values):
    html = "<table>\n"
    if len(x_labels) > 0:
        html += "<tr>\n"
        for label in x_labels:
            html += "\t<th>%s</th>" % label
        html += "</tr>\n"
    for i, label in enumerate(y_labels):
        if i % 2 == 0:
            html += "<tr class='even'>\n\t<td>%s</td>\n" % label
            for j in range(len(values)):
                html += "\t<td class='even'>%s</td>\n" % values[j][i]
        else:
            html += "<tr>\n\t<td>%s</td>\n" % label
            for j in range(len(values)):
                html += "\t<td>%s</td>\n" % values[j][i]
        html += "</tr>\n"
    html += "</table>\n"
    return html

def load_insert_distribution(data):
    html = """
    <h4>Insert Size Distribution</h4>
    <div class="chart" id="insert"></div>
    <script>
    google.load('visualization', '1', {packages: ['corechart', 'line']});
    google.setOnLoadCallback(drawBackgroundColor);

    function drawBackgroundColor() {
          var data = new google.visualization.DataTable();
          data.addColumn('number', 'X');
          data.addColumn('number', 'estimate');
          data.addColumn('number', 'observed');

          data.addRows([
    """
    insert_dist = data['insert_distribution'][...]
    Y = insert_dist[:, 0]
    X = insert_dist[:, 1].astype(numpy.float64)
    temp = numpy.log(X[1:])
    X_span = numpy.exp(numpy.mean(temp[1:]-temp[:-1]) / 2.0)
    X[-1] = X[-2] * X_span
    X[0] = X[1] / X_span
    X[1:-1] *= X_span
    best_insert, est_Y = learn_GMM(numpy.log(X[1:]), Y[1:])
    est_Y = numpy.r_[0., est_Y]
    for i in range(X.shape[0]):
        html += "[%f, %f, %i]," % (X[i], est_Y[i], Y[i])
        if i % 10 == 9:
            html += "\n\t\t"
    html += """
    ]);
    var options = {
        hAxis: {
            title: 'Insert size',
            logScale: true
        },
        vAxis: {
            title: 'Read count'
        },
        legend: {
            position: 'none'
        },
    };
    var chart = new google.visualization.LineChart(document.getElementById('insert'));
    chart.draw(data, options);
    }
    </script>
    """
    return best_insert, html

def load_cis_interaction_distribution(data):
    html = """
    <h4>Fend Intra-chromosomal Interaction Distribution</h4>
    <div class="chart" id="cis"></div>
    <script>
    google.load('visualization', '1', {packages: ['corechart', 'line']});
    google.setOnLoadCallback(drawBackgroundColor);

    function drawBackgroundColor() {
          var data = new google.visualization.DataTable();
          data.addColumn('number', 'X');
          data.addColumn('number', 'Y');

          data.addRows([
    """
    Y = data['cis_interaction_distribution'][...]
    X = numpy.arange(Y.shape[0]).astype(numpy.float64)
    X[0] = 0.1
    for i in range(X.shape[0]):
        html += "[%f, %i]," % (X[i], Y[i])
        if i % 6 == 5:
            html += "\n\t\t"
    html += """
    ]);
    var options = {
        hAxis: {
            title: 'Number of intra-chromosomal interaction partners',
            logScale: true
        },
        vAxis: {
            title: 'Number of fends'
        },
        legend: {
            position: 'none'
        },
    };
    var chart = new google.visualization.LineChart(document.getElementById('cis'));
    chart.draw(data, options);
    }
    </script>
    """
    return html

def load_trans_interaction_distribution(data):
    html = """
    <h4>Fend Inter-chromosomal Interaction Distribution</h4>
    <div class="chart" id="trans"></div>
    <script>
    google.load('visualization', '1', {packages: ['corechart', 'line']});
    google.setOnLoadCallback(drawBackgroundColor);

    function drawBackgroundColor() {
          var data = new google.visualization.DataTable();
          data.addColumn('number', 'X');
          data.addColumn('number', 'Y');

          data.addRows([
    """
    Y = data['trans_interaction_distribution'][...]
    X = numpy.arange(Y.shape[0]).astype(numpy.float64)
    X[0] = 0.1
    for i in range(X.shape[0]):
        html += "[%f, %i]," % (X[i], Y[i])
        if i % 6 == 5:
            html += "\n\t\t"
    html += """
    ]);
    var options = {
        hAxis: {
            title: 'Number of inter-chromosomal interaction partners',
            logScale: true
        },
        vAxis: {
            title: 'Number of fends'
        },
        legend: {
            position: 'none'
        },
    };
    var chart = new google.visualization.LineChart(document.getElementById('trans'));
    chart.draw(data, options);
    }
    </script>
    """
    return html

def int2str(i):
    string = str(i)
    j = (len(string) - 1) % 3 + 1
    formatted = string[:j]
    while j < len(string):
        formatted += ',' + string[j:(j + 3)]
        j += 3
    return formatted

def learn_GMM(X, Y):
    try:
        num_dists = 6
        temp = numpy.copy(Y).astype(numpy.float64)
        for i in range(1, Y.shape[0]):
            temp[i] += temp[i - 1]
        temp /= temp[-1]
        indices = numpy.searchsorted(temp, numpy.arange(num_dists + 1) / float(num_dists + 1))
        mu = numpy.zeros(num_dists, dtype=numpy.float64)
        sigma = numpy.zeros(num_dists, dtype=numpy.float64)
        M = numpy.zeros(num_dists, dtype=numpy.float64)
        for i in range(num_dists):
            mu[i] = (numpy.sum(Y[indices[i]:indices[i + 1]] * X[indices[i]:indices[i + 1]]) /
                     numpy.sum(Y[indices[i]:indices[i + 1]]))
            sigma[i] = (X[indices[i + 1]] - X[indices[i]]) / 4.0
            M[i] = 1.0 / num_dists
        prob = numpy.zeros((X.shape[0], M.shape[0]), dtype=numpy.float64)
        for i in range(1000):
            valid = numpy.where(M > 0.01)[0]
            prob.fill(0.0)
            for j in valid:
                #prob[:, j] = numpy.exp(-(X - mu[j]) ** 2 / (2 * sigma[j])) / ((sigma[j] * 2 * numpy.pi) ** 0.5) * M[j]
                prob[:, j] = numpy.exp(-0.5 * (numpy.log(2.0 * numpy.pi * sigma[j]) + (X - mu[j]) ** 2.0 / sigma[j])) * M[j]
            prob /= numpy.sum(prob, axis=1).reshape(-1, 1)
            prob *= Y.reshape(-1, 1)
            M = numpy.sum(prob, axis=0)
            mu[valid] = numpy.sum(prob[:, valid] * X.reshape(-1, 1), axis=0) / M[valid]
            sigma[valid] = numpy.sum(prob[:, valid] * (X.reshape(-1, 1) - mu[valid]) ** 2, axis=0) / M[valid]
            M /= numpy.sum(M)
        est_Y = (numpy.exp(-(X.reshape(-1, 1) - mu.reshape(1, -1)) ** 2 / (2 * sigma.reshape(1, -1))) /
                          ((sigma.reshape(1, -1) * 2 * numpy.pi) ** 0.5) * M.reshape(1, -1))
        est_Y = numpy.sum(est_Y, axis=1)
        est_Y *= numpy.sum(Y) / numpy.sum(est_Y)
        peaks = numpy.where(est_Y[1:-1] > numpy.maximum(est_Y[:-2], est_Y[2:]))[0] + 1
        troughs = numpy.r_[0, numpy.where(est_Y[1:-1] < numpy.minimum(est_Y[:-2], est_Y[2:]))[0] + 1,
                           est_Y.shape[0] - 1]
        indices = numpy.searchsorted(troughs, peaks)
        peak_scores = est_Y[peaks] - numpy.maximum(est_Y[troughs[indices]], est_Y[troughs[indices - 1]])
        where1 = numpy.where(peak_scores == numpy.amax(peak_scores))[0][0]
        peak1 = peaks[where1]
        peak_scores[where1] = -numpy.inf
        if peak_scores.shape[0] > 1:
            where2 = numpy.where(peak_scores == numpy.amax(peak_scores))[0][0]
            peak2 = peaks[where2]
            while max(X[peak1], X[peak2]) < numpy.log(300) and numpy.amax(peak_scores) > -numpy.inf:
                if peak1 < peak2:
                    peak1, peak2 = peak2, peak1
                    where1, where2 = where2, where1
                peak_scores[where2] = -numpy.inf
                where2 = numpy.where(peak_scores == numpy.amax(peak_scores))[0][0]
                peak2 = peaks[where2]
            if numpy.amax(peak_scores) > -numpy.inf:
                if peak2 < peak1:
                    peak1, peak2 = peak2, peak1
                split = numpy.exp(X[numpy.where(est_Y[peak1:peak2] == numpy.amin(est_Y[peak1:peak2]))[0][0] + peak1])
            else:
                split = 650
        else:
            split = 650
        return split, est_Y
    except:
        return 650, numpy.zeros(Y.shape[0], dtype=numpy.float32)

if __name__ == "__main__":
    main()
