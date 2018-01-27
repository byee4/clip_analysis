"""
Just a helper class that helps visualize 
"""
from matplotlib.colors import rgb2hex, colorConverter
from collections import defaultdict

class Clusters(dict):
    """
    Converts rgb values to hex for viewing.
    In reality I don't need this and can just use a dict, but
    useful for viewing stuff.
    """

    def _repr_html_(self):
        html = '<table style="border: 0;">'
        for c in self:
            hx = rgb2hex(colorConverter.to_rgb(c))
            html += '<tr style="border: 0;">' \
                    '<td style="background-color: {0}; ' \
                    'border: 0;">' \
                    '<code style="background-color: {0};">'.format(hx)
            html += c + '</code></td>'
            html += '<td style="border: 0"><code>'
            html += repr(self[c]) + '</code>'
            html += '</td></tr>'

        html += '</table>'

        return html


def get_cluster_classes(den, label='ivl'):
    """
    Appends gene names to their corresponding cluster ids
    """
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = Clusters()
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes