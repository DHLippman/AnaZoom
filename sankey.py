"""
Author:         David Henry Lippman
File:           main.py
Date created:   03/20/20
Date modified:  03/20/20

"""

import numpy as np
import plotly.graph_objects as go


trials = 1e3
foPass = np.array([3, 63, 112, 78, 1])
foFail = trials - sum(foPass)
rtPass = np.array([0, 7, 23, 1, 0])
rtFail = foPass - rtPass

# node = dict(pad=25,
#             thickness=30,
#             line=dict(color=None),
#             label=["16 million Monte Carlo trials",     # 0
#                    "Invalid first-order solutions",     # 1
#                    "PPNP first-order solutions (3)",    # 2
#                    "PNPP first-order solutions (63)",   # 3
#                    "PNPN first-order solutions (112)",  # 4
#                    "PNNP first-order solutions (78)",   # 5
#                    "NPPN first-order solutions (1)",    # 6
#                    "Invalid ray traceable solutions",   # 7
#                    "PNPP ray traceable solutions (7)",  # 8
#                    "PNPN ray traceable solutions (23)", # 9
#                    "PNNP ray traceable solutions (1)",  # 10
#                    ],
#             color=["blue", "red", "green", "green", "green", "green", "green",
#                    "red", "green", "green", "green"])

node = dict(pad=25,
            thickness=30,
            line=dict(color=None),
            label=["16 million Monte Carlo trials",     # 0
                   "Invalid solutions",     # 1
                   "PPNP solutions (3)",    # 2
                   "PNPP solutions (63)",   # 3
                   "PNPN solutions (112)",  # 4
                   "PNNP solutions (78)",   # 5
                   "NPPN solutions (1)",    # 6
                   "Invalid solutions",     # 7
                   "PNPP solutions (7)",    # 8
                   "PNPN solutions (23)",   # 9
                   "PNNP solutions (1)",    # 10
                   ],
            color=["blue", "red", "green", "green", "green", "green", "green",
                   "red", "green", "green", "green"])

link = dict(source=[0, 0, 0, 0, 0, 0,
                    1, 2, 3, 3, 4, 4, 5, 5, 6],
            target=[1, 2, 3, 4, 5, 6,
                    7, 7, 7, 8, 7, 9, 7, 10, 7],
            value=[foFail, foPass[0], foPass[1], foPass[2], foPass[3], foPass[4],
                   foFail, rtFail[0], rtFail[1], rtPass[1], rtFail[2], rtPass[2], rtFail[3], rtPass[3], rtFail[4]])

data=[go.Sankey(node=node, link=link)]

fig = go.Figure(data=data)

fig.update_layout(title_text="Monte Carlo First-Order Search",
                  font=dict(family="Arial", size=30))
fig.show()