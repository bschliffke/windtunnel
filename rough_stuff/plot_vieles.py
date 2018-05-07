def plot_convergence_test(data1,data2,data3,data4,data5,data6,data7,data8,data9):
    """ Plots results of convergence tests performed on each quantity
    into one 3x3 plot. This is a very limited function and is only intended to
    give a brief overview of the convergence test results using dictionaries as
    input objects.
    @parameter: data[0,...,9], type = dictionary"""
    
    fig, ax = plt.subplots(3,3, figsize=(12,7))
    for i, key in enumerate([key for key in data1.keys()]):
        ax[0,0].plot([i] * len(data1[key]), data1[key], color='navy',
                      linestyle='None',marker='o', markersize=15)                  
    del i, key
    for i, key in enumerate([key for key in data2.keys()]):
        ax[0,1].plot([i] * len(data2[key]), data2[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data3.keys()]):
        ax[0,2].plot([i] * len(data3[key]), data3[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data4.keys()]):
        ax[1,0].plot([i] * len(data4[key]), data4[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data5.keys()]):
        ax[1,1].plot([i] * len(data5[key]), data5[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data6.keys()]):
        ax[1,2].plot([i] * len(data6[key]), data6[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data7.keys()]):
        ax[2,0].plot([i] * len(data7[key]), data7[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data8.keys()]):
        ax[2,1].plot([i] * len(data8[key]), data8[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    for i, key in enumerate([key for key in data9.keys()]):
        ax[2,2].plot([i] * len(data9[key]), data9[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
    del i, key
    
    for i in range(ax.shape[0]):
        for k in range(ax.shape[1]):
            ax[i][k].tick_params(labelsize=14)
        #    ax[k].spines['left'].set_position('zero')
            ax[i][k].spines['right'].set_color('none') # .set_visible(False)
        #    ax[k].spines['bottom'].set_position('zero')
            ax[i][k].spines['top'].set_color('none') # .set_visible(False)
            ax[i][k].xaxis.set_ticks_position('bottom')  
            ax[i][k].yaxis.set_ticks_position('left')
    del i,k 
    
    ax[0,0].set(xticks=np.arange(len(data1.keys())+1),
              xticklabels=[key for key in data1.keys()],
              xlim=(-0.5, len(data1.keys())-0.5))
    ax[0,1].set(xticks=np.arange(len(data2.keys())+1),
              xticklabels=[key for key in data2.keys()],
              xlim=(-0.5, len(data2.keys())-0.5))
    ax[0,2].set(xticks=np.arange(len(data3.keys())+1),
              xticklabels=[key for key in data3.keys()],
              xlim=(-0.5, len(data3.keys())-0.5))
    ax[1,0].set(xticks=np.arange(len(data4.keys())+1),
              xticklabels=[key for key in data4.keys()],
              xlim=(-0.5, len(data4.keys())-0.5))
    ax[1,1].set(xticks=np.arange(len(data5.keys())+1),
              xticklabels=[key for key in data5.keys()],
              xlim=(-0.5, len(data5.keys())-0.5))
    ax[1,2].set(xticks=np.arange(len(data6.keys())+1),
              xticklabels=[key for key in data6.keys()],
              xlim=(-0.5, len(data6.keys())-0.5))
    ax[2,0].set(xticks=np.arange(len(data7.keys())+1),
              xticklabels=[key for key in data7.keys()],
              xlim=(-0.5, len(data7.keys())-0.5))
    ax[2,1].set(xticks=np.arange(len(data8.keys())+1),
              xticklabels=[key for key in data8.keys()],
              xlim=(-0.5, len(data8.keys())-0.5))
    ax[2,2].set(xticks=np.arange(len(data9.keys())+1),
              xticklabels=[key for key in data9.keys()],
              xlim=(-0.5, len(data9.keys())-0.5))