# The following code changes default matplotlib.pyplot parameters to standardize plots between notebooks
def change_defaults():
    import matplotlib.pyplot as plt_now
    plt_now.rcParams['font.family'] = ['Arial']
    plt_now.rcParams['text.color'] = 'xkcd:dark gray'
    plt_now.rcParams['axes.edgecolor'] = 'xkcd:dark gray'
    plt_now.rcParams['axes.labelcolor'] = 'xkcd:dark gray'
    plt_now.rcParams['axes.spines.top'] = False
    plt_now.rcParams['axes.spines.right'] = False
    plt_now.rcParams['xtick.labelcolor'] = 'xkcd:dark gray'
    plt_now.rcParams['ytick.labelcolor'] = 'xkcd:dark gray'
    plt_now.rcParams['legend.frameon'] = False
    plt_now.rcParams['legend.fontsize'] = 14
    plt_now.rcParams['xtick.labelsize'] = 18
    plt_now.rcParams['ytick.labelsize'] = 18
    plt_now.rcParams['scatter.edgecolors'] = 'xkcd:dark gray'
    plt_now.rcParams['savefig.format'] = 'pdf'
    plt_now.rcParams['savefig.transparent'] = True
    plt_now.rcParams['legend.handletextpad'] = 0
    plt_now.rcParams['xtick.color'] = 'xkcd:dark gray'
    plt_now.rcParams['ytick.color'] = 'xkcd:dark gray'