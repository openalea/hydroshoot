#!/usr/bin/env python2
# -*- coding: utf-8 -*-


def plot_figure_6():
    """Generates figure 6 of the paper. This figure traces the hourly values of absorbed irradiance, leaf temeprature,
    net plant carbon assimilation rate, and net plant transpiration rate.
    """

    training_color = {'gdc': 'blue', 'vsp': 'red', 'lyre': 'green'}

    beg_date = datetime(2009, 7, 29, 00, 00, 0, )
    end_date = datetime(2009, 7, 29, 23, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='H')

    meteo_df = pd.read_csv(example_pth / 'virtual_canopies' / 'gdc' / 'meteo.input',
                           sep=';', decimal='.', index_col='time')  # all simus have the same meteo data
    meteo_df.index = pd.DatetimeIndex(meteo_df.index)
    meteo_df = meteo_df.loc[datet]

    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharex=True, figsize=(6.69, 6))
    [ax.grid() for ax in axs.flatten()]

    axs[1, 0].plot(datet, meteo_df['Tac'], 'k--')

    for training in ('gdc', 'vsp', 'lyre'):
        pth = example_pth / 'virtual_canopies' / training

        sims_df = pd.read_csv(pth / 'output' / 'time_series.output', sep=';', decimal='.', index_col='time')
        sims_df.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims_df.index]

        axs[0, 0].plot(datet, sims_df['Rg'], label=training, color=training_color[training])
        axs[0, 1].plot(datet, sims_df['An'], label=training, color=training_color[training])
        axs[1, 0].plot(datet, sims_df['Tleaf'], label=training, color=training_color[training])
        axs[1, 1].plot(datet, sims_df['E'], label=training, color=training_color[training])

    # some layout
    for iax, ax in enumerate(axs.flatten()):
        ax.text(0.875, 0.9, ('(a)', '(b)', '(c)', '(d)')[iax],
                transform=ax.transAxes, fontdict={'size': 11})

    axs[0, 0].set(xlim=(beg_date, end_date), ylim=(0, 600),
                  ylabel='$\mathregular{\Phi_{R_g, plant}\/[W\/m^{-2}_{ground}]}$')
    axs[1, 0].set(xlabel='hour', ylim=(10, 40),
                  ylabel='$\mathregular{Temperature\/[^\circ C]}$')
    axs[0, 1].set(ylim=(-10, 35),
                  ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$')
    axs[1, 1].set(xlabel='hour', ylim=(0, 400),
                  ylabel='$\mathregular{E_{plant}\/[g\/h^{-1}]}$')

    for ax in axs[1, :]:
        ax.set(xlim=(beg_date, end_date))
        ax.set_xticklabels(datet, rotation=90)
        ax.xaxis.set_major_locator(dates.HourLocator())
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))

    axs[0, 0].legend(loc='upper left', prop={'size': 11})
    handles, _ = axs[1, 0].get_legend_handles_labels()
    axs[1, 0].legend(handles=[handles[0], ], labels=('$\mathregular{T_{air}}$',), loc='upper left',
                     prop={'size': 11})

    axs[1, 0].set_xticks(axs[1, 0].get_xticks()[:-1:2])

    fig.tight_layout()

    fig.savefig('fig_6.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_7():
    """Generates figure 7 of the paper. This figure shows a snapshot at solar midday of water potential distribution
    across three virtual canopies.
    """
    training_color = {'gdc': 'blue', 'vsp': 'red', 'lyre': 'green'}

    obs_date = datetime(2009, 7, 29, 14, 00, 0, )

    fig, axs = pyplot.subplots(nrows=3, ncols=2, sharex=True, figsize=(6.69, 6),
                               gridspec_kw={'width_ratios': [0.8, 0.2]})

    for i, training in enumerate(('vsp', 'gdc', 'lyre')):
        pth = example_pth / 'virtual_canopies' / training

        g, _ = mtg_load(str(pth) + '/output/', 'mtg%s' % obs_date.strftime('%Y%m%d%H%M%S'))

        axs[i, 0] = display.property_map(g, 'psi_head', color=training_color[training],
                                         ax=axs[i, 0])
        axs[i, 0].set(ylim=(0, 250), xlim=(-1.7, -1.2))
        axs[i, 0].legend([training], prop={'size': 11})
        axs[i, 0].xaxis.labelpad = 5

        axs[i, 1].boxplot(g.property('psi_head').values(), vert=False, sym='')

    # some layout
    [ax.set_xlabel('') for ax in axs[:2, 0]]

    for ax in axs[:, 1]:
        ax.tick_params(axis='y', which='both', labelleft='off')
        ax.set_xticklabels(np.arange(-1.7, -1.1, 0.1), rotation=90)

    axs[2, 1].set_xlabel('$\mathregular{\Psi\/[MPa]}$')

    fig.tight_layout()

    fig.savefig('fig_7.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_8():
    """Generates figure 8 of the paper. This figure shows a snapshot at solar midday of leaf temperature distribution
    across three virtual canopies.
    """
    training_color = {'gdc': 'blue', 'vsp': 'red', 'lyre': 'green'}

    obs_date = datetime(2009, 7, 29, 14, 00, 0, )

    fig, axs = pyplot.subplots(nrows=3, ncols=2, sharex=True, figsize=(6.69, 6),
                               gridspec_kw={'width_ratios': [0.8, 0.2]})

    for i, training in enumerate(('vsp', 'gdc', 'lyre')):
        pth = example_pth / 'virtual_canopies' / training

        g, _ = mtg_load(str(pth) + '/output/', 'mtg%s' % obs_date.strftime('%Y%m%d%H%M%S'))

        axs[i, 0] = display.property_map(g, 'Tlc', color=training_color[training],
                                         ax=axs[i, 0])
        axs[i, 0].set(ylim=(0, 250))
        axs[i, 0].legend([training], prop={'size': 11})
        axs[i, 0].xaxis.labelpad = 5

        axs[i, 1].boxplot(g.property('Tlc').values(), vert=False, sym='')

    # some layout
    [ax.set_xlabel('') for ax in axs[:2, 0]]

    for ax in axs[:, 1]:
        ax.tick_params(axis='y', which='both', labelleft='off')
        ax.set_xticklabels(range(32, 42), rotation=90)

    axs[2, 1].set_xlabel('$\mathregular{T_{leaf}\/[^{\circ}C]}$')

    fig.tight_layout()

    fig.savefig('fig_8.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_9():
    """Generates figure 9 of the paper. This figure compares, simulated to observed xylem water potential.
    """

    beg_date = datetime(2012, 8, 1, 00, 00, 0, )
    end_date = datetime(2012, 8, 3, 00, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='D')

    fig, axs = pyplot.subplots(nrows=3, ncols=3, sharex='all', sharey='all', figsize=(6.69, 6))

    for it, training in enumerate(('gdc_can1_grapevine', 'gdc_can2_grapevine', 'gdc_can3_grapevine')):

        pth = example_pth / training

        for i_day, date in enumerate(datet):
            ax = axs[it, i_day]
            obs_date = date + pd.Timedelta(hours=13)
            g, _ = mtg_load(str(pth) + '/output/', 'mtg%s' % obs_date.strftime('%Y%m%d%H%M%S'))
            ax = display.property_map(g, 'psi_head', ax=ax, prop2='Eabs', color='grey',
                                      colormap='autumn', add_color_bar=False)

            obs_df = pd.read_csv(pth / 'var.obs', sep=';', index_col='date')
            obs_df.index = pd.DatetimeIndex(obs_df.index)
            psi_leaf = obs_df.loc[date, 'psi_leaf']
            ax.add_patch(patches.Rectangle((max(psi_leaf), 50),
                                           min(psi_leaf) - max(psi_leaf), 250,
                                           color='0.8', zorder=-1))
            if i_day == 0:
                ax.text(0.05, 0.075, 'Canopy%d' % (it + 1),
                        transform=ax.transAxes, fontdict={'size': 11})

            if it == 0:
                ax.text(0.675, 0.825, obs_date.strftime('%-d %b'),
                        transform=ax.transAxes, fontdict={'size': 11})

    [axi.set_xticklabels(axi.get_xticks(), rotation=90) for axi in axs[2]]

    for ax, axi in enumerate(axs):
        if ax < 2:
            [axi[i_day].set_xlabel('') for i_day in range(3)]
        [axi[i_day].set_ylabel('') for i_day in (1, 2)]
        [axi[i_day].legend_.remove() for i_day in range(3)]

    fig.tight_layout()

    fig.subplots_adjust(bottom=0.3)
    cbar_ax = fig.add_axes([0.395, 0.1, 0.30, 0.04])
    norm = colors.Normalize(0, vmax=2000)
    cbar = colorbar.ColorbarBase(cbar_ax, cmap='autumn', orientation='horizontal', norm=norm)
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    cbar.set_label('$\mathregular{PPFD\/[\mu mol\/m^{-2}_{leaf}\/s^{-1}]}$', labelpad=-20, x=-0.4)

    fig.savefig('fig_9.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_10():
    """Generates figure 10 of the paper. This figure compares, simulated to observed stomatal conductance rates.
    """

    beg_date = datetime(2012, 8, 1, 00, 00, 0, )
    end_date = datetime(2012, 8, 3, 00, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='D')

    fig, axs = pyplot.subplots(nrows=3, ncols=3, sharex='all', sharey='all', figsize=(6.69, 6))

    for it, training in enumerate(('gdc_can1_grapevine', 'gdc_can2_grapevine', 'gdc_can3_grapevine')):

        pth = example_pth / training

        for i_day, date in enumerate(datet):
            ax = axs[it, i_day]
            obs_date = date + pd.Timedelta(hours=13)
            g, _ = mtg_load(str(pth) + '/output/', 'mtg%s' % obs_date.strftime('%Y%m%d%H%M%S'))
            ax = display.property_map(g, 'gs', ax=ax, prop2='Eabs', color='grey',
                                      colormap='autumn', add_color_bar=False)

            obs_df = pd.read_csv(pth / 'var.obs', sep=';', index_col='date')
            obs_df.index = pd.DatetimeIndex(obs_df.index)
            gs_leaf = obs_df.loc[date, 'gs'] / 1000.
            ax.add_patch(patches.Rectangle((max(gs_leaf), 50),
                                           min(gs_leaf) - max(gs_leaf), 250,
                                           color='0.8', zorder=-1))
            if i_day == 0:
                ax.text(0.05, 0.075, 'Canopy%d' % (it + 1),
                        transform=ax.transAxes, fontdict={'size': 11})

            if it == 0:
                ax.text(0.675, 0.825, obs_date.strftime('%-d %b'),
                        transform=ax.transAxes, fontdict={'size': 11})

    [axi.set_xticklabels(axi.get_xticks(), rotation=90) for axi in axs[2]]

    for ax, axi in enumerate(axs):
        if ax < 2:
            [axi[i_day].set_xlabel('') for i_day in range(3)]
        [axi[i_day].set_ylabel('') for i_day in (1, 2)]
        [axi[i_day].legend_.remove() for i_day in range(3)]

    fig.tight_layout()

    fig.subplots_adjust(bottom=0.3)
    cbar_ax = fig.add_axes([0.395, 0.1, 0.30, 0.04])
    norm = colors.Normalize(0, vmax=2000)
    cbar = colorbar.ColorbarBase(cbar_ax, cmap='autumn', orientation='horizontal', norm=norm)
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=90)
    cbar.set_label('$\mathregular{PPFD\/[\mu mol\/m^{-2}_{leaf}\/s^{-1}]}$', labelpad=-20, x=-0.4)

    fig.savefig('fig_10.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_11():
    """Generates figure 11 of the paper. This figure compares simulated to observed leaf temperatures.
    """
    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharey='row', figsize=(6.69, 6))

    for iax, training in enumerate(('vsp_ww_grapevine', 'vsp_ws_grapevine')):
        pth = example_pth / training
        axs[0, iax].grid()
        axs[1, iax].grid()

        # set dates ************************
        beg_date = datetime(2009, 7, 29, 00, 00, 0, )
        end_date = datetime(2009, 8, 1, 23, 00, 0, )
        datet = pd.date_range(beg_date, end_date, freq='H')

        # read observations
        obs = pd.read_csv(pth / 'temp.obs', sep=';', decimal='.', index_col='date')
        obs.index = [datetime.strptime(s, "%d/%m/%Y %H:%M") for s in obs.index]

        # read simulations
        sims = pd.read_csv(pth / 'output' / 'time_series.output',
                           sep=';', decimal='.', index_col='time')
        sims.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims.index]
        sims.index = pd.DatetimeIndex(sims.index)

        # plot median simulated leaf temperature
        axs[0, iax].plot(sims['Tleaf'], '-k', label='$\mathregular{T_{leaf}}$', linewidth=1)

        # plot simulated temperature magnitude
        med_sim = np.array([])
        q1_sim = np.array([])
        q3_sim = np.array([])
        for date in datet:
            g, _ = mtg_load(str(pth) + '/output/', 'mtg%s' % date.strftime('%Y%m%d%H%M%S'))
            leaf_temp_sim = g.property('Tlc').values()
            q1_sim = np.append(q1_sim, min(leaf_temp_sim))
            q3_sim = np.append(q3_sim, max(leaf_temp_sim))
            med_sim = np.append(med_sim, np.median(leaf_temp_sim))
            print(date)
        axs[0, iax].fill_between(datet, q1_sim, q3_sim, color='red', alpha=0.5, zorder=0)

        # plot observed temperature magnitude
        med_obs = np.array([])
        q1_obs = np.array([])
        q3_obs = np.array([])
        for row, date in enumerate(datet):
            pos = date.toordinal() + date.hour / 24.
            leaf_temp_obs = obs.loc[date, ['Tleaf%d' % d for d in (2, 3, 4, 5, 6, 8, 9, 10)]]
            leaf_temp_obs = leaf_temp_obs[~np.isnan(leaf_temp_obs)]
            axs[0, iax].boxplot(leaf_temp_obs, positions=[pos], widths=0.01)
            q1_obs = np.append(q1_obs, min(leaf_temp_obs))
            q3_obs = np.append(q3_obs, max(leaf_temp_obs))
            med_obs = np.append(med_obs, np.median(leaf_temp_obs))
            print(date)

        # compare obs to sim temperature on 1:1 plot
        axs[1, iax].plot((10, 50), (10, 50), 'k--')
        axs[1, iax].errorbar(x=med_obs, y=med_sim,
                             xerr=(np.nan_to_num(med_obs - q1_obs), np.nan_to_num(q3_obs - med_obs)),
                             yerr=(np.nan_to_num(med_sim - q1_sim), np.nan_to_num(q3_sim - med_sim)),
                             fmt='ro', ecolor='0.5', capthick=1)

        # plot MBE and RMSE results on it
        diff_obs_sim = (med_sim - med_obs)[~np.isnan(med_sim - med_obs)]
        bias = diff_obs_sim.mean()
        rmse = np.sqrt(diff_obs_sim ** 2).mean()

        axs[1, iax].text(0.50, 0.20, '$\mathregular{MBE\/=\/%.3f}$' % bias,
                         transform=axs[1, iax].transAxes, fontdict={'size': 11})
        axs[1, iax].text(0.50, 0.10, '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                         transform=axs[1, iax].transAxes, fontdict={'size': 11})

    # some layout

    axs[0, 0].set_ylabel('$\mathregular{Temperature\/[^{\circ}C]}$')
    axs[1, 0].set_ylabel('$\mathregular{T_{leaf, sim}\/[^{\circ}C]}$')

    for ax_upper, ax_lower in zip(axs[0, :], axs[1, :]):
        ax_upper.set(xlim=(beg_date, end_date), ylim=(5, 50), xlabel='')
        ax_upper.set_xticklabels(datet, rotation=45)
        ax_upper.xaxis.set_major_locator(dates.DayLocator())
        ax_upper.xaxis.set_major_formatter(dates.DateFormatter('%-d %b'))
        ax_lower.set(xlim=(10, 50), ylim=(10, 50),
                     xlabel='$\mathregular{T_{leaf, obs}\/[^{\circ}C]}$')

    for iax, ax in enumerate(axs.flatten()):
        ax.text(0.05, 0.875, ('(a)', '(c)', '(b)', '(d)')[iax],
                transform=ax.transAxes, fontdict={'size': 11})

    fig.tight_layout()
    fig.savefig('fig_11.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_12():
    """Generates figure 12 of the paper. This figure compares, leaf-to-leaf, simulated to observed leaf and leaf-to-air
    temperature.
    """

    fig, axs = pyplot.subplots(ncols=3, figsize=(6.69, 2.7))
    [ax.grid() for ax in axs]

    daily_temp_obs = np.array([])
    daily_temp_sim = np.array([])
    daily_temp_air = np.array([])

    for iax, training in enumerate(('vsp_ww_grapevine', 'vsp_ws_grapevine')):

        pth = example_pth / training

        # set dates ************************
        beg_date = datetime(2009, 7, 29, 00, 00, 0, )
        end_date = datetime(2009, 8, 1, 23, 00, 0, )
        datet = pd.date_range(beg_date, end_date, freq='H')

        # read observations
        obs = pd.read_csv(pth / 'temp.obs', sep=';', decimal='.', index_col='date')
        obs.index = [datetime.strptime(s, "%d/%m/%Y %H:%M") for s in obs.index]

        # read simulations
        sims = pd.read_csv(pth / 'output' / 'time_series.output',
                           sep=';', decimal='.', index_col='time')
        sims.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims.index]
        sims.index = pd.DatetimeIndex(sims.index)

        # plot simulated temperature magnitude
        med_sim = np.array([])
        q_1_sim = np.array([])
        q_3_sim = np.array([])
        for date in datet:
            g, _ = mtg_load(str(pth) + '/output/', 'mtg%s' % date.strftime('%Y%m%d%H%M%S'))
            leaf_temp_sim = g.property('Tlc').values()
            q_1_sim = np.append(q_1_sim, min(leaf_temp_sim))
            q_3_sim = np.append(q_3_sim, max(leaf_temp_sim))
            med_sim = np.append(med_sim, np.median(leaf_temp_sim))
            print(date)

        # get observed temperature magnitude
        med_obs = np.array([])
        q1_obs = np.array([])
        q3_obs = np.array([])
        for row, date in enumerate(datet):
            pos = date.toordinal() + date.hour / 24.
            leaf_temp_obs = obs.loc[date, ['Tleaf%d' % d for d in (2, 3, 4, 5, 6, 8, 9, 10)]]
            leaf_temp_obs = leaf_temp_obs[~np.isnan(leaf_temp_obs)]
            q1_obs = np.append(q1_obs, min(leaf_temp_obs))
            q3_obs = np.append(q3_obs, max(leaf_temp_obs))
            med_obs = np.append(med_obs, np.median(leaf_temp_obs))
            print(date)

        # read meteo data
        meteo_df = pd.read_csv(pth / 'meteo.input', sep=';', decimal='.', index_col='time')
        meteo_df.index = pd.DatetimeIndex(meteo_df.index)

        air_temperature = meteo_df.loc[datet, 'Tac']

        daily_temp_obs = np.append(daily_temp_obs, q1_obs)
        daily_temp_obs = np.append(daily_temp_obs, q3_obs)

        daily_temp_sim = np.append(daily_temp_sim, q_1_sim)
        daily_temp_sim = np.append(daily_temp_sim, q_3_sim)

        daily_temp_air = np.append(daily_temp_air, air_temperature)
        daily_temp_air = np.append(daily_temp_air, air_temperature)

        # minimum temperature
        axs[0].plot(q1_obs - air_temperature, q_1_sim - air_temperature,
                    ['b^', 'r^'][iax], label=['WW', 'WD'][iax])
        # maximum temperature
        axs[0].plot(q3_obs - air_temperature, q_3_sim - air_temperature,
                    ['bo', 'ro'][iax], label=['WW', 'WD'][iax])

        axs[1].plot(q1_obs, q_1_sim, ['b^', 'r^'][iax], label=['WW', 'WD'][iax])
        axs[1].plot(q3_obs, q_3_sim, ['bo', 'ro'][iax], label=['WW', 'WD'][iax])

        axs[2].plot(q3_obs - q1_obs, q_3_sim - q_1_sim, ['bs', 'rs'][iax], label=['WW', 'WD'][iax])

    # some layout
    axs[0].plot((-8, 12), (-8, 12), 'k--')
    axs[0].set(xlabel='$\mathregular{T_{leaf, obs}-T_{air}\/[^\circ C]}$',
               ylabel='$\mathregular{T_{leaf, sim}-T_{air}\/[^\circ C]}$',
               xlim=(-8, 12), ylim=(-8, 12))

    axs[1].plot((10, 45), (10, 45), 'k--')
    axs[1].set(xlabel='$\mathregular{T_{leaf, obs}\/[^\circ C]}$',
               ylabel='$\mathregular{T_{leaf, sim}\/[^\circ C]}$')
    axs[1].xaxis.set_major_locator(ticker.MultipleLocator(20))

    axs[2].plot((-2, 14), (-2, 14), 'k--')
    axs[2].set(xlabel='$\mathregular{\Delta T_{leaf, obs}\/[^\circ C]}$',
               ylabel='$\mathregular{\Delta T_{leaf, sim}\/[^\circ C]}$')
    axs[2].xaxis.set_major_locator(ticker.IndexLocator(base=2, offset=2))
    axs[2].xaxis.set_major_locator(ticker.MultipleLocator(5))

    for i in range(3):
        if i == 0:
            x, y = daily_temp_obs - daily_temp_air, daily_temp_sim - daily_temp_air
        elif i == 1:
            x, y = daily_temp_obs, daily_temp_sim
        else:
            x = np.array([])
            x = np.append(x, daily_temp_obs[96:2 * 96] - daily_temp_obs[:96])
            x = np.append(x, daily_temp_obs[3 * 96:4 * 96] - daily_temp_obs[2 * 96:3 * 96])
            y = np.array([])
            y = np.append(y, daily_temp_sim[96:2 * 96] - daily_temp_sim[:96])
            y = np.append(y, daily_temp_sim[3 * 96:4 * 96] - daily_temp_sim[2 * 96:3 * 96])

        diff_y_x = (y - x)[~np.isnan(y - x)]
        bias = (diff_y_x).mean()
        rmse = np.sqrt((diff_y_x ** 2).mean())

        axs[i].text(0.05, -0.6, '$\mathregular{MBE\/=\/%.3f}$' % bias,
                    transform=axs[i].transAxes, fontdict={'size': 11})
        axs[i].text(0.05, -0.7, '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                    transform=axs[i].transAxes, fontdict={'size': 11})

        axs[i].text(0.05, 0.875, ['(a)', '(b)', '(c)'][i],
                    transform=axs[i].transAxes, fontdict={'size': 11})

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.4)
    fig.savefig('fig_12.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_13():
    """Generates figure 13 of the paper.
    """
    # set dates
    beg_date = datetime(2009, 7, 29, 00, 00, 0, )
    end_date = datetime(2009, 8, 1, 23, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='H')

    fig = pyplot.figure()
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(left=0.135, right=0.65, top=0.975, bottom=0.125, wspace=0.1, hspace=0.35)
    ax1 = pyplot.subplot(gs1[0, 0])
    ax2 = pyplot.subplot(gs1[0, 1], sharex=ax1)
    ax3 = pyplot.subplot(gs1[1, 0], sharex=ax1)
    ax4 = pyplot.subplot(gs1[1, 1], sharex=ax1)

    gs2 = gridspec.GridSpec(2, 1)
    gs2.update(left=0.67, right=0.865, top=0.975, bottom=0.125, hspace=0.35)
    ax5 = pyplot.subplot(gs2[0])
    ax6 = pyplot.subplot(gs2[1])

    fig.set_figheight(6)
    fig.set_figwidth(6.69)
    axs = np.array([[ax1, ax2, ax5],
                    [ax3, ax4, ax6]])

    [ax.grid() for ax in axs.flatten()]

    for iax, training in enumerate(('vsp_ww_grapevine', 'vsp_ws_grapevine')):
        pth = example_pth / training

        # read observations
        obs_tab = pd.read_csv(pth / 'gas.obs', sep=';', decimal='.', header=0)
        obs_tab.time = pd.DatetimeIndex(obs_tab.time)
        obs_tab = obs_tab.set_index(obs_tab.time)

        axs[0, iax].plot(obs_tab['time'], obs_tab['An_plante'],
                         color='0.8', marker='o', markeredgecolor='none', label='$\mathregular{A_{n,\/obs}}$')
        axs[1, iax].plot(obs_tab['time'], obs_tab['E_plante'],
                         color='0.8', marker='o', markeredgecolor='none', label='$\mathregular{E_{obs}}$')

        # read simulations
        sims = pd.read_csv(pth / 'output' / 'time_series.output',
                           sep=';', decimal='.', index_col='time')
        sims.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims.index]
        sims.index = pd.DatetimeIndex(sims.index)

        axs[0, iax].plot(sims['An'],
                         'b-', label='$\mathregular{A_{n,\/sim}}$', linewidth=1)
        axs[1, iax].plot(sims['E'],
                         'b-', label='$\mathregular{E_{sim}}$', linewidth=1)

        sims['itime'] = sims.index
        obs_tab['time'] = pd.DatetimeIndex(obs_tab.time)

        time_group = []
        for itime in obs_tab.time:
            if itime.minute < 30:
                time = itime - pd.Timedelta(minutes=itime.minute)
            else:
                time = itime + pd.Timedelta(minutes=60 - itime.minute)
            time_group.append(time)
        obs_tab['itime'] = time_group

        obs_grouped = obs_tab.groupby('itime').aggregate(np.mean)
        obs_grouped['itime'] = obs_grouped.index

        m_df = pd.merge(obs_grouped, sims)

        axs[0, 2].plot(m_df['An_plante'], m_df['An'], ('bo', 'ro')[iax])
        axs[1, 2].plot(m_df['E_plante'], m_df['E'], ('bo', 'ro')[iax])

        for egas, igas in enumerate(('An', 'E')):
            x = m_df[igas + '_plante'].values
            y = m_df[igas].values

            xindex = np.isfinite(x)
            yindex = np.isfinite(y)
            xyindex = xindex * yindex
            x, y = x[xyindex], y[xyindex]

            bias = (y - x).mean()
            rmse = np.sqrt(((x - y) ** 2).mean())

            axs[egas, iax].text(0.25, 0.90,
                                '$\mathregular{MBE\/=\/%.1f}$' % bias,
                                transform=axs[egas, iax].transAxes,
                                fontdict={'size': 11})
            axs[egas, iax].text(0.25, 0.8,
                                '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                                transform=axs[egas, iax].transAxes,
                                fontdict={'size': 11})

    axs[0, 2].plot((-20, 50), (-20, 50), 'k--')
    axs[1, 2].plot((-200, 1000), (-200, 1000), 'k--')

    for i, ax in enumerate(axs.flatten()):
        ax.text(0.05, 0.9, ('(a)', '(b)', '(e)', '(c)', '(d)', '(f)')[i],
                transform=ax.transAxes, fontdict={'size': 11})

    # some layout
    axs[0, 0].set(
        ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$',
        xlim=(beg_date, end_date),
        ylim=(-20, 50))
    axs[1, 0].set(
        ylabel='$\mathregular{E_{plant}\/[g\/h^{-1}]}$',
        xlim=(beg_date, end_date),
        ylim=(-200, 1000))
    axs[0, 1].set(
        xlim=(beg_date, end_date),
        ylim=(-20, 50))
    axs[1, 1].set(
        xlim=(beg_date, end_date),
        ylim=(-200, 1000))
    axs[0, 2].set(
        xlabel='$\mathregular{A_{n, plant, obs}\/[\mu mol\/s^{-1}]}$',
        ylabel='$\mathregular{A_{n, plant, sim}\/[\mu mol\/s^{-1}]}$',
        xlim=(-20, 50),
        ylim=(-20, 50))
    axs[1, 2].set(
        xlabel='$\mathregular{E_{plant, obs}\/[g\/h^{-1}]}$',
        ylabel='$\mathregular{E_{plant, sim}\/[g\/h^{-1}]}$',
        xlim=(-200, 1000),
        ylim=(-200, 1000))

    axs[0, 1].get_yaxis().set_ticklabels('')
    axs[1, 1].get_yaxis().set_ticklabels('')
    axs[1, 1].get_yaxis().set_ticklabels('')
    axs[0, 2].yaxis.tick_right()
    axs[1, 2].yaxis.tick_right()
    axs[0, 2].yaxis.set_label_position('right')
    axs[1, 2].yaxis.set_label_position('right')
    axs[0, 2].xaxis.set_major_locator(ticker.MultipleLocator(25))
    axs[1, 2].xaxis.set_major_locator(ticker.MultipleLocator(500))

    for ax in axs[:, :2].flatten():
        ax.set_xticklabels(datet, rotation=90)
        ax.xaxis.set_major_locator(dates.DayLocator())
        ax.xaxis.set_major_formatter(dates.DateFormatter('%-d %b'))

    fig.savefig('fig_13.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_14():
    """Generates figure 14 of the paper. This figure compares, simulated to observed plant transpiration rates.
    """

    fig, axs = pyplot.subplots(nrows=3, ncols=4, sharey='row', figsize=(6.69, 6))

    for i_treat, training in enumerate(('gdc_can1_grapevine', 'gdc_can2_grapevine', 'gdc_can3_grapevine')):

        pth = example_pth / training

        # read observations
        obs_df = pd.read_csv(pth / 'sapflow.obs', sep=';', decimal='.', index_col='date')
        obs_df.index = [datetime.strptime(s, "%d/%m/%Y %H:%M") for s in obs_df.index]

        time_group = []
        for itime in obs_df.index:
            if itime.minute < 30:
                time = itime - pd.Timedelta(minutes=itime.minute)
            else:
                time = itime + pd.Timedelta(minutes=60 - itime.minute)
            time_group.append(time)
        obs_df['itime'] = time_group

        obs_grouped = obs_df.groupby('itime').aggregate(np.mean)
        obs_grouped['itime'] = obs_grouped.index

        # read simulations
        sims = pd.read_csv(pth / 'output' / 'time_series.output',
                           sep=';', decimal='.', index_col='time')
        sims.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims.index]
        sims.index = pd.DatetimeIndex(sims.index)

        sims['itime'] = sims.index
        m_df = pd.merge(obs_grouped, sims)

        x = m_df['west'].values + m_df['east'].values
        y = m_df['arm2'].values + m_df['arm1'].values

        x_index = np.isfinite(x)
        y_index = np.isfinite(y)
        xy_index = x_index * y_index
        x, y = x[xy_index], y[xy_index]

        for iax, ax in enumerate(axs[i_treat, :]):
            ax.xaxis.grid(which='minor', zorder=0)
            ax.yaxis.grid(which='major', zorder=0)

            ax.plot(sims.index, sims['arm1'], c='r', linewidth=1)
            ax.plot(sims.index, sims['arm2'], c='b', linewidth=1)

            ax.plot(obs_df.index, obs_df['east'], label='East', color='red',
                    marker='o', markeredgecolor='none', alpha=0.1)
            ax.plot(obs_df.index, obs_df['west'], label='West', color='blue',
                    marker='o', markeredgecolor='none', alpha=0.1)

            x_t = x[iax * 24:iax * 24 + 23]
            y_t = y[iax * 24:iax * 24 + 23]

            x_index = np.isfinite(x_t)
            y_index = np.isfinite(y_t)
            xy_index = x_index * y_index
            x_t, y_t = x_t[xy_index], y_t[xy_index]

            bias = (y_t - x_t).mean()
            rmse = np.sqrt(((x_t - y_t) ** 2).mean())

            ax.text(0.45, 0.675,
                    '$\mathregular{MBE\/=\/%.1f}$' % bias,
                    transform=ax.transAxes, fontdict={'size': 7})
            ax.text(0.45, 0.575,
                    '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                    transform=ax.transAxes, fontdict={'size': 7})

    # some layout
    for day in range(4):
        day_sdate = datetime(2012, 8, 1, 00, 00, 0, ) + timedelta(days=day)
        day_edate = day_sdate + timedelta(hours=23)
        for can in range(3):
            axs[can, day].set_xlim(day_sdate, day_edate)
            axs[can, day].yaxis.set_major_locator(ticker.MultipleLocator(100))
            if can != 2:
                axs[can, day].xaxis.set_ticklabels('')
                axs[can, day].xaxis.set_minor_locator(dates.HourLocator(interval=5))
            else:
                axs[can, day].set_xticklabels(pd.date_range(day_sdate, day_edate, freq='H'),
                                              rotation=90)
                axs[can, day].xaxis.set_major_locator(dates.DayLocator())
                axs[can, day].xaxis.set_major_formatter(dates.DateFormatter('%-d %b'))
                axs[can, day].xaxis.set_minor_locator(dates.HourLocator(interval=5))
                axs[can, day].xaxis.set_minor_formatter(dates.DateFormatter('%H'))
                axs[can, day].tick_params(axis='x', which='major', pad=15)

    axs[1, 0].set_ylabel('$\mathregular{E\/[g\/h^{-1}]}$')

    axs[0, 0].set_ylim(0, 800)
    axs[1, 0].set_ylim(0, 450)
    axs[2, 0].set_ylim(0, 500)

    for can in range(3):
        axs[can, 0].text(0.05, 0.85, 'Canopy%s' % str(can + 1),
                         transform=axs[can, 0].transAxes, fontdict={'size': 11})

    fig.tight_layout()
    fig.subplots_adjust(left=0.13, bottom=0.3)

    h1, l1 = axs[2, 2].get_legend_handles_labels()

    axs[2, 0].legend(h1, ('$\mathregular{SapEast_{sim}}$', '$\mathregular{SapWest_{sim}}$',
                          '$\mathregular{SapEast_{obs}}$', '$\mathregular{SapWest_{obs}}$'),
                     frameon=True, bbox_to_anchor=(-0.36, -1.2, 2, .102),
                     loc=3, ncol=8, prop={'size': 11})
    fig.subplots_adjust(wspace=0.075)

    fig.savefig('fig_14.png', dpi=600.)
    pyplot.close(fig)


def plot_figure_15():
    """Generates figure 15 of the paper. This figure compares simulated to observed plant photosynthesis and
    transpiration rates using 4 simulation parameters' sets (modularity test)
    """
    # set dates
    beg_date = datetime(2009, 7, 29, 00, 00, 0, )
    end_date = datetime(2009, 8, 1, 23, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='H')

    style = ('b-', 'k--', 'k-.', 'k-', 'k:')
    vpd_air = np.vectorize(VPDa)

    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharey='row', sharex='all', figsize=(6.69, 6))
    [ax.grid() for ax in axs.flatten()]

    for i_treat, treat in enumerate(('ww', 'ws')):

        pth = example_pth / ('vsp_%s_grapevine' % treat)

        # read and plot observations
        obs_df = pd.read_csv(pth / 'gas.obs', sep=';', decimal='.')
        obs_df.time = pd.DatetimeIndex(obs_df.time)
        obs_df = obs_df.set_index(obs_df.time)

        axs[0, i_treat].plot(obs_df['time'], obs_df['An_plante'],
                             color='0.8', marker='o', markeredgecolor='none', label='$\mathregular{A_{n,\/obs}}$')
        axs[1, i_treat].plot(obs_df['time'], obs_df['E_plante'],
                             color='0.8', marker='o', markeredgecolor='none', label='$\mathregular{E_{obs}}$')

        # read and plot vpd
        meteo_df = pd.read_csv(pth / 'meteo.input', sep=';', decimal='.', index_col='time')
        meteo_df.index = pd.DatetimeIndex(meteo_df.index)
        meteo_df['vpd'] = vpd_air(meteo_df.Tac, meteo_df.Tac, meteo_df.hs)
        axt = axs[1, i_treat].twinx()
        axt.plot(datet, meteo_df.loc[datet, 'vpd'], 'r-')

        # read and plot simulations
        for case in range(5):
            if case == 0:
                sim_pth = pth / 'output' / 'time_series.output'
            else:
                sim_pth = example_pth / 'modularity' / treat / ('sim_%d' % case) / 'output' / 'time_series.output'

            sims = pd.read_csv(sim_pth, sep=';', decimal='.', index_col='time')
            sims.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims.index]
            sims.index = pd.DatetimeIndex(sims.index)

            axs[0, i_treat].plot(sims['An'], style[case], label='sim0', linewidth=1)
            axs[1, i_treat].plot(sims['E'], style[case], label='sim0', linewidth=1)

    for i, ax in enumerate(axs.flatten()):
        ax.text(0.05, 0.9, ('(a)', '(b)', '(c)', '(d)')[i], transform=ax.transAxes,
                fontdict={'size': 11})

    axs[0, 0].set(
        ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$',
        ylim=(-20, 50))
    axs[1, 0].set(
        ylabel='$\mathregular{E_{plant}\/[g\/h^{-1}]}$',
        xlim=(beg_date, end_date),
        ylim=(-200, 1600))

    [ax.set_ylim(0, 5) for ax in fig.get_axes()[-2:]]
    fig.get_axes()[-1].set(ylabel='VPD [kPa]')

    for ax in axs[1, :].flatten():
        ax.set_xticklabels(datet, rotation=90)
        ax.xaxis.set_major_locator(dates.DayLocator())
        ax.xaxis.set_major_formatter(dates.DateFormatter('%-d %b'))

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.25)

    h1, l1 = axs[1, 1].get_legend_handles_labels()
    h2, l2 = fig.get_axes()[-1].get_legend_handles_labels()

    axs[1, 0].legend(h1 + h2, ['obs'] + ['sim%d' % d for d in range(5)] + ['VPD'], frameon=True,
                     bbox_to_anchor=(0.15, -0.75, 2, .102), loc=3, ncol=4,
                     prop={'size': 11})

    fig.savefig('fig_15.png', dpi=600.)
    pyplot.close(fig)


def write_table_1():
    """Generates table 1 of the paper. This table compares simulated to observed plant photosynthesis and
    transpiration rates using 4 simulation parameters' sets (modularity test)
    """
    indices = [['ww'] * 5 + ['ws'] * 5,
               ['sim%d' % d for d in range(5)] * 2]
    multiple_index = pd.MultiIndex.from_tuples(list(zip(*indices)), names=['treat', 'sim'])

    df = pd.DataFrame(columns=('an_mbe', 'an_rmse', 'e_mbe', 'e_rmse'),
                      index=multiple_index)

    for i_treat, treat in enumerate(('ww', 'ws')):

        pth = example_pth / ('vsp_%s_grapevine' % treat)

        # read observations
        obs_df = pd.read_csv(pth / 'gas.obs', sep=';', decimal='.')
        obs_df.time = pd.DatetimeIndex(obs_df.time)
        obs_df = obs_df.set_index(obs_df.time)

        # read simulations
        for case in range(5):
            if case == 0:
                sim_pth = pth / 'output' / 'time_series.output'
            else:
                sim_pth = example_pth / 'modularity' / treat / ('sim_%d' % case) / 'output' / 'time_series.output'

            sims_df = pd.read_csv(sim_pth, sep=';', decimal='.', index_col='time')
            sims_df.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims_df.index]
            sims_df.index = pd.DatetimeIndex(sims_df.index)

            # match sims to obs
            sims_df['itime'] = sims_df.index
            obs_df['time'] = pd.DatetimeIndex(obs_df.time)

            time_group = []
            for itime in obs_df.index:
                if itime.minute < 30:
                    time = itime - pd.Timedelta(minutes=itime.minute)
                else:
                    time = itime + pd.Timedelta(minutes=60 - itime.minute)
                time_group.append(time)
            obs_df['itime'] = time_group

            obs_grouped = obs_df.groupby('itime').aggregate(np.mean)
            obs_grouped['itime'] = obs_grouped.index

            m_df = pd.merge(obs_grouped, sims_df)

            for gas in ('An', 'E'):
                rate_obs = m_df['%s_plante' % gas].values
                rate_sim = m_df[gas].values

                x_index = np.isfinite(rate_obs)
                y_index = np.isfinite(rate_sim)
                xy_index = x_index * y_index
                rate_obs, rate_sim = rate_obs[xy_index], rate_sim[xy_index]

                bias = (rate_sim - rate_obs).mean()
                rmse = np.sqrt(((rate_obs - rate_sim) ** 2).mean())

                df.loc[(treat, 'sim%d' % case), '%s_mbe' % gas.lower()] = bias
                df.loc[(treat, 'sim%d' % case), '%s_rmse' % gas.lower()] = rmse

    print(df)
    df.to_csv('table_1.csv')


def estimate_energy_balance_contribution():
    """Generates table 1 of the paper. This table compares simulated to observed plant photosynthesis and
    transpiration rates using 4 simulation parameters' sets (modularity test)
    """
    indices = [['ww'] * 4 + ['ws'] * 4,
               ['day%d' % d for d in range(1, 5)] * 2]
    multiple_index = pd.MultiIndex.from_tuples(list(zip(*indices)), names=['treat', 'day'])

    df = pd.DataFrame(columns=('e_obs', 'e_sim0', 'e_sim3'),
                      index=multiple_index)

    for i_treat, treat in enumerate(('ww', 'ws')):

        pth = example_pth / ('vsp_%s_grapevine' % treat)

        # read observations
        obs_df = pd.read_csv(pth / 'gas.obs', sep=';', decimal='.')
        obs_df.time = pd.DatetimeIndex(obs_df.time)
        obs_df = obs_df.set_index(obs_df.time)

        # read simulations
        for case in (0, 3):
            if case == 0:
                sim_pth = pth / 'output' / 'time_series.output'
            else:
                sim_pth = example_pth / 'modularity' / treat / ('sim_%d' % case) / 'output' / 'time_series.output'

            sims_df = pd.read_csv(sim_pth, sep=';', decimal='.', index_col='time')
            sims_df.index = [datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in sims_df.index]
            sims_df.index = pd.DatetimeIndex(sims_df.index)

            # match sims to obs
            sims_df['itime'] = sims_df.index
            obs_df['time'] = pd.DatetimeIndex(obs_df.time)

            time_group = []
            for itime in obs_df.index:
                if itime.minute < 30:
                    time = itime - pd.Timedelta(minutes=itime.minute)
                else:
                    time = itime + pd.Timedelta(minutes=60 - itime.minute)
                time_group.append(time)
            obs_df['itime'] = time_group

            obs_grouped = obs_df.groupby('itime').aggregate(np.mean)
            obs_grouped['itime'] = obs_grouped.index

            m_df = pd.merge(obs_grouped, sims_df)

            rate_obs = m_df['E_plante'].values
            rate_sim = m_df['E'].values

            for day in range(4):
                rate_obs_daily = rate_obs[day * 24: day * 24 + 24]
                rate_sim_daily = rate_sim[day * 24: day * 24 + 24]

                x_index = np.isfinite(rate_obs_daily)
                y_index = np.isfinite(rate_sim_daily)
                xy_index = x_index * y_index
                rate_obs_daily, rate_sim_daily = rate_obs_daily[xy_index], rate_sim_daily[xy_index]

                df.loc[(treat, 'day%d' % (day + 1)), 'e_obs'] = sum(rate_obs_daily)
                df.loc[(treat, 'day%d' % (day + 1)), 'e_sim%d' % case] = sum(rate_sim_daily)

    df['energy_effect'] = (df['e_sim3'] - df['e_sim0']) / df['e_sim0']

    print(df)


if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    from pathlib import Path
    from datetime import datetime, timedelta
    from matplotlib import dates, pyplot, patches, colors, colorbar, rcParams, ticker, gridspec

    from hydroshoot.architecture import mtg_load
    from hydroshoot.utilities import vapor_pressure_deficit as VPDa
    from hydroshoot import display

    rcParams.update({'font.size': 11})
    pyplot.style.use('seaborn-ticks')

    example_pth = Path(__file__).parents[2] / 'example'

    plot_figure_6()
    plot_figure_7()
    plot_figure_8()
    plot_figure_9()
    plot_figure_10()
    plot_figure_11()
    plot_figure_12()
    plot_figure_13()
    plot_figure_14()
    plot_figure_15()
    write_table_1()
    estimate_energy_balance_contribution()
