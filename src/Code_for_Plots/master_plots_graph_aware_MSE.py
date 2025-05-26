'''
Plot MSE
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
    x_label = ['ratio', 'tp', 'size', 'percent']
    x_var = ['ratio', 'p', 'n', 'pct']
    x_plot = ['$r$', '$p$', '$n$','$pct$']
    x_label = ['tp', 'size']
    x_var = ['p', 'n']
    x_plot = ['$p$', '$n$']

    graph = "srgg"
    for beta in [2]:
        title = ['$\\beta='+str(beta)+'$','$\\beta='+str(beta)+', n=10000, r=2$','$\\beta='+str(beta)+', p=0.2, r=2$']
        #title = ['$\\beta='+str(beta)+', n=10000, p=0.2, r=2$']
        est_names = ['Reg', 'VIM', 'SNIPE('+str(beta)+')', 'Lin\'s', 'DM']
        for ind in [0,1,2]:
            plot(graph,x_var[ind],x_label[ind],'deg'+str(beta),x_plot[ind],title[ind],est_names,permute=True)


def plot(graph,x_var,x_label,model,x_plot,title,est_names,permute=False):
    experiment = '-'+x_label+'-'+model
    print(experiment)

    # Create and save plots
    df = pd.read_csv(load_path+graph+experiment+'-SNIPE.csv')

    # Compute MSE
    df["biassq"] = df["Bias_squared"]
    # print(df)
    df2 = df.groupby([x_var,'Estimator']).agg('mean', numeric_only=True)
    print(df2["biassq"])

    # Uncomment line below for LaTeX font
    plt.rc('text', usetex=True)
    
    # Plot with all the estimators
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # '#2077b4' - blue; '#ff7f0e' - orange; '#2ba02c' - green; '#d62728' - red; '#9468bd' - purple
    sns.lineplot(x=x_var, y='biassq', hue='Estimator', style='Estimator', data=df2, legend='brief', markers=True, style_order=est_names, hue_order= est_names, palette=['#ff7f0e', '#2ba02c', '#9468bd', '#d62728', '#2077b4', '#FFC0CB'])

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
    # Add zoomed-in inset
    axins = inset_axes(ax, width="95%", height="95%", loc='upper left', 
                   bbox_to_anchor=(0, 0.5, .5, .5), bbox_transform=ax.transAxes)

    sns.lineplot(
        x=x_var, y='biassq', hue='Estimator', style='Estimator',
        data=df2.reset_index(),
        hue_order=est_names, style_order=est_names,
        palette=['#ff7f0e', '#2ba02c', '#9468bd', '#d62728', '#2077b4', '#FFC0CB'], markers=False, ax=axins, legend=False
    )

    # Adjust zoom region here
    if x_var == 'n': 
        axins.set_xlim(8000, 10010)
        axins.set_ylim(10, 60)
    elif x_var == 'p':
        axins.set_xlim(0.099, 0.155)
        axins.set_ylim(50, 120)
    axins.set_xticks([])
    axins.set_yticks([])
    axins.set_xlabel("")
    axins.set_ylabel("")

    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="gray")
    
    if model == 'deg1':
        pass
        ax.set_ylim(0,50)
    else:
        pass
        #ax.set_ylim(0,100)
    ax.set_xlabel(x_plot, fontsize = 18)
    ax.set_ylabel("MSE", fontsize = 18)
    ax.set_title(title, fontsize=20)
    handles, labels = ax.get_legend_handles_labels()

    if permute:
        order = [4,0,1,3,2] #1 - DM75, 0- DM, 3-LSprop, 4-SNIpe, 2-LSnum
        ax.legend(handles=[handles[i] for i in order], labels=[labels[i] for i in order], loc='upper right', fontsize = 14)
    else:
        ax.legend(handles=handles, labels=labels, loc='upper right', fontsize = 14)

    plt.savefig(save_path+graph+experiment+'_SNIPE_MSE.pdf')
    plt.close()

if __name__ == "__main__":
    main()
