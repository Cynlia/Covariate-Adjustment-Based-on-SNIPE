'''
Plot relative bias and experimental standard deviation
'''

# Setup
from matplotlib import rcParams
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Run with these settings to reproduce the MSE plots from the origianl paper
load_path = 'outputFiles/graph_aware/'
save_path = 'outputFiles/graph_aware/'

# Run with these settings to produce MSE plots from new data
#load_path = 'outputFiles/new/'
#save_path = 'outputFiles/new/'

def main():    
    x_label = ['ratio', 'tp', 'size']
    x_var = ['ratio', 'p', 'n']
    x_plot = ['$r$', '$p$', '$n$']
    x_label = ['tp', 'size']
    x_var = ['p', 'n']
    x_plot = ['$p$', '$n$']
    #x_label = ['percent']
    #x_var = ['pct']
    #x_plot = ['$pct$']
    graph = "srgg"
    for beta in [2]:
        if graph == "sw":
            title = ['$\\beta='+str(beta)+', n=9216, p=0.2$','$\\beta='+str(beta)+', n=9216, r=2$','$\\beta='+str(beta)+', p=0.2, r=2$']
        else:
            title = ['$\\beta='+str(beta)+'$','$\\beta='+str(beta)+'$','$\\beta='+str(beta)+', p=0.2, r=2$']
            #title = ['$\\beta='+str(beta)+', n=10000, p=0.2, r=2$']
        est_names = ['Reg', 'VIM', 'SNIPE('+str(beta)+')', 'Lin\'s', 'DM']
        for ind in [0,1,2]:
            plot(graph,x_var[ind],x_label[ind],'deg'+str(beta),x_plot[ind],title[ind],est_names,permute=True)


def plot(graph,x_var,x_label,model,x_plot,title,est_names,permute=False):
    experiment = '-'+x_label+'-'+model
    print(experiment)

    # Create and save plots
    df = pd.read_csv(load_path+graph+experiment+'-SNIPE.csv')

    plt.rc('text', usetex=True)
    
    # Plot with all the estimators
    fig = plt.figure()
    ax = fig.add_subplot(111)

    sns.lineplot(x=x_var, y='Relative_Bias', hue='Estimator', style='Estimator', data=df, ci='sd', legend='brief', markers=True, style_order=est_names, hue_order= est_names, palette=['#ff7f0e', '#2ba02c', '#9468bd', '#d62728', '#2077b4', '#FFC0CB'])
    ax.set_ylim(-1,1)
    ax.set_xlabel(x_plot, fontsize = 18)
    ax.set_ylabel("Relative Bias", fontsize = 18)
    ax.set_title(title, fontsize=20)
    handles, labels = ax.get_legend_handles_labels()

    if permute:
        order = [4,0,1,3,2]
        ax.legend(handles=[handles[i] for i in order], labels=[labels[i] for i in order], loc='upper right', fontsize = 14)
    else:
        ax.legend(handles=handles, labels=labels, loc='upper right', fontsize = 14)

    plt.savefig(save_path+graph+experiment+'_SNIPE.pdf')
    plt.close()

if __name__ == "__main__":
    main()
