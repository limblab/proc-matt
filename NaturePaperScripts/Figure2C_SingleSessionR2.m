clear
clc;
close all;

dataSummary;


% the good potent/null sessions
sessions = { ...
    'Chewie','2016-10-07'; ...
    };

basenames = {'trainad','trainad','trainad'};
extranames = {'','',''};
array_models = {'M1-M1','PMd-PMd','PMd-M1'};

pert            = 'FF';
tasks           = {'CO'};
dates           = sessions(:,2);
monkeys         = unique(sessions(:,1));

which_metric    = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff      = 0.01;
pr2_op          = 'min'; % which operation for filtering ('min','max','mean','median')
pr2_ad_check    = false; % only keeps cells that predict in WO
do_same_cells   = false; % really only works for testing tweaks of model

plot_op         = 'mean';
group_size      = 20;
how_to_group    = 'slide'; % average, pool, slide

do_norm         = true;
do_diff         = true;
remove_outliers = true;
num_outlier_std = 10;

error_bars      = 'ste'; % 'boot','ste'
num_bootstraps  = 1000;

do_reg_line     = false;
add_plot        = true;
do_subplot      = false;
filter_trials   = false;

epochs          = {'AD'};%{'BL','AD','WO'};

if do_norm
    min_y = -1.3; max_y = 0.1;
else
    min_y = -0.1; max_y = 0.1;
end

glm_r2_plots;
