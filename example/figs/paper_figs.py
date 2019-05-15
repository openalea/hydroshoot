#!/usr/bin/env python2
# -*- coding: utf-8 -*-


def plot_figure_13():
    """Generates figure 13 of the paper.
    """
    # set dates
    beg_date = datetime(2009, 7, 29, 00, 00, 0, )
    end_date = datetime(2009, 8, 1, 23, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='H')

    fig, axs = pyplot.subplots(nrows=2, ncols=3, figsize=(10, 6))
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

            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            r2 = round(r_value ** 2, 3)
            bias = (y - x).mean()
            rmse = np.sqrt(((x - y) ** 2).mean())

            axs[egas, iax].text(0.65, 0.90,
                                '$\mathregular{MBE\/=\/%.1f}$' % bias,
                                transform=axs[egas, iax].transAxes,
                                fontdict={'size': 7})
            axs[egas, iax].text(0.65, 0.8,
                                '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                                transform=axs[egas, iax].transAxes,
                                fontdict={'size': 7})

    axs[0, 2].plot((-20, 50), (-20, 50), 'k--')
    axs[1, 2].plot((-200, 1000), (-200, 1000), 'k--')

    for i, ax in enumerate(axs.flatten()):
        ax.text(0.05, 0.9, ('(a)', '(b)', '(e)', '(c)', '(d)', '(f)')[i],
                transform=ax.transAxes)

    axs[0, 0].set(
        xlabel='Date',
        ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$',
        xlim=(beg_date, end_date),
        ylim=(-20, 50))
    axs[1, 0].set(
        xlabel='Date',
        ylabel='$\mathregular{E_{plant}\/[g\/h^{-1}]}$',
        xlim=(beg_date, end_date),
        ylim=(-200, 1000))
    axs[0, 1].set(
        xlabel='Date',
        ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$',
        xlim=(beg_date, end_date), ylim=(-20, 50))
    axs[1, 1].set(
        xlabel='Date',
        ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$',
        xlim=(beg_date, end_date), ylim=(-200, 1000))
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

    for ax in axs[:, :2].flatten():
        ax.set_xticklabels(datet, rotation=90)
        ax.xaxis.set_major_locator(dates.DayLocator())
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d %m'))

    fig.tight_layout()
    fig.savefig('fig_13.png')

    pyplot.show(fig)


def plot_meteo_vsp():
    beg_date = datetime(2009, 7, 29, 00, 00, 0, )
    end_date = datetime(2009, 8, 1, 23, 00, 0, )
    datet = pd.date_range(beg_date, end_date, freq='H')

    fig, axs = pyplot.subplots(nrows=3)
    [ax.grid() for ax in axs.flatten()]

    for iax, training in enumerate(('vsp_ww_grapevine', 'vsp_ws_grapevine')):
        pth = example_pth / training

        # read meteo
        obs_tab = pd.read_csv(pth / 'meteo.input', sep=';', decimal='.', header=0)
        obs_tab.time = pd.DatetimeIndex(obs_tab.time)
        obs_tab = obs_tab.set_index(obs_tab.time)
        obs_tab = obs_tab.loc[datet]
        obs_tab['depression'] = ([0] * 10 + [1] * 5 + [0] * 9) * 4

        axs[0].plot(obs_tab['time'], obs_tab['Tac'],
                    ('b-', 'r-')[iax], label=('ww', 'wd')[iax])
        axs[1].plot(obs_tab[obs_tab['depression'] == 1]['Rg'], obs_tab[obs_tab['depression'] == 1]['Tac'],
                    ('bo', 'ro')[iax])
    axs[2].plot(obs_tab['time'], obs_tab['Rg'])
    pyplot.show(fig)


def plot_figure_11():
    """Generates figure 11 of the paper. This figure compares simulated to observed leaf temperatures.
    """
    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharey='row')

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
            g, _ = mtg_load(str(pth) + '/output/mtg', date.strftime('%Y%m%d%H%M%S'))
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

        axs[1, iax].text(0.05, 0.80, '$\mathregular{MBE\/=\/%.3f}$' % bias,
                         transform=axs[1, iax].transAxes, fontdict={'size': 7})
        axs[1, iax].text(0.05, 0.75, '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                         transform=axs[1, iax].transAxes, fontdict={'size': 7})

    # some layout

    axs[0, 0].set_ylabel('$\mathregular{Temperature\/[^{\circ}C]}$')
    axs[1, 0].set_ylabel('$\mathregular{T_{leaf, sim}\/[^{\circ}C]}$')

    for ax_upper, ax_lower in zip(axs[0, :], axs[1, :]):
        ax_upper.set(xlim=(beg_date, end_date), ylim=(5, 50), xlabel='Date')
        ax_upper.set_xticklabels(datet, rotation=45)
        ax_upper.xaxis.set_major_locator(dates.DayLocator())
        ax_upper.xaxis.set_major_formatter(dates.DateFormatter('%d %m'))
        ax_lower.set(xlim=(10, 50), ylim=(10, 50),
                     xlabel='$\mathregular{T_{leaf, obs}\/[^{\circ}C]}$')

    for iax, ax in enumerate(axs.flatten()):
        ax.text(0.05, 0.9, ('(a)', '(c)', '(b)', '(d)')[iax],
                transform=ax.transAxes, fontdict={'size': 10})

    fig.tight_layout()
    fig.savefig('fig_11.png')

    pyplot.show(fig)


def plot_figure_12():
    """Generates figure 12 of the paper. This figure compares, leaf-to-leaf, simulated to observed leaf and leaf-to-air
    temperature.
    """

    fig, axs = pyplot.subplots(ncols=3, figsize=(13, 4.5))
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
            g, _ = mtg_load(str(pth) + '/output/mtg', date.strftime('%Y%m%d%H%M%S'))
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

    axs[2].plot((-2, 14), (-2, 14), 'k--')
    axs[2].set(xlabel='$\mathregular{\Delta T_{leaf, obs}\/[^\circ C]}$',
               ylabel='$\mathregular{\Delta T_{leaf, sim}\/[^\circ C]}$')

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

        axs[i].text(0.05, 0.80, '$\mathregular{MBE\/=\/%.3f}$' % bias,
                    transform=axs[i].transAxes, fontdict={'size': 7})
        axs[i].text(0.05, 0.75, '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                    transform=axs[i].transAxes, fontdict={'size': 7})

        axs[i].text(0.05, 0.9, ['(a)', '(b)', '(c)'][i],
                    transform=axs[i].transAxes, fontdict={'size': 10})

    fig.tight_layout()
    fig.savefig('fig_12.png')

    pyplot.show(fig)


def plot_figure_14():
    """Generates figure 14 of the paper. This figure compares, simulated to observed plant transpiration rates.
    """

    fig, axs = pyplot.subplots(nrows=3, ncols=4, sharex='col', sharey='row', figsize=(10, 6))

    for i_treat, training in enumerate(('gdc_can1_grapevine', 'gdc_can2_grapevine', 'gdc_can3_grapevine')):

        pth = example_pth / training

        beg_date = datetime(2012, 8, 01, 00, 00, 0, )
        end_date = datetime(2012, 8, 04, 23, 00, 0, )
        datet = pd.date_range(beg_date, end_date, freq='H')

        # read observations
        obs_df = pd.read_csv(pth / 'sapflow.obs', sep=';', decimal='.', index_col='date')
        obs_df.index = [datetime.strptime(s, "%d/%m/%Y %H:%M") for s in obs_df.index]
        # obs_df.index = pd.DatetimeIndex(obs_df.idnex)

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
        y = m_df['sapWest'].values + m_df['sapEast'].values

        x_index = np.isfinite(x)
        y_index = np.isfinite(y)
        xy_index = x_index * y_index
        x, y = x[xy_index], y[xy_index]

        for iax, ax in enumerate(axs[i_treat, :]):
            ax.xaxis.grid(which='minor', zorder=0)
            ax.yaxis.grid(which='major', zorder=0)

            ax.plot(sims.index, sims.sapEast, c='r', linewidth=1)
            ax.plot(sims.index, sims.sapWest, c='b', linewidth=1)

            ax.plot(obs_df.index, obs_df['east'], label='East', color='red',
                    marker='o', markeredgecolor='none', alpha=0.1)
            ax.plot(obs_df.index, obs_df['west'], label='West', color='blue',
                    marker='o', markeredgecolor='none', alpha=0.1)

            x_t = x[i_treat * 24:i_treat * 24 + 23]
            y_t = y[i_treat * 24:i_treat * 24 + 23]

            x_index = np.isfinite(x_t)
            y_index = np.isfinite(y_t)
            xy_index = x_index * y_index
            x_t, y_t = x_t[xy_index], y_t[xy_index]

            bias = (y_t - x_t).mean()
            rmse = np.sqrt(((x_t - y_t) ** 2).mean())

            ax.text(0.45, 0.875,
                    '$\mathregular{MBE\/=\/%.1f}$' % bias,
                    transform=ax.transAxes, fontdict={'size': 7})
            ax.text(0.45, 0.775,
                    '$\mathregular{RMSE\/=\/%.1f}$' % rmse,
                    transform=ax.transAxes, fontdict={'size': 7})

    # some layout
    for day in range(4):
        day_sdate = datetime(2012, 8, 01, 00, 00, 0, ) + timedelta(days=day)
        day_edate = day_sdate + timedelta(hours=23)
        for can in range(3):
            axs[can, day].set_xlim(day_sdate, day_edate)

    axs[1, 0].set_ylabel('$\mathregular{E\/[g\/h^{-1}]}$')

    axs[0, 0].set_ylim(0, 800)
    axs[1, 0].set_ylim(0, 450)
    axs[2, 0].set_ylim(0, 500)

    for can in range(3):
        axs[can, 3].text(0.5, 0.5, 'Canopy%s' % str(can + 1),
                         transform=axs[can, 3].transAxes, fontdict={'size': 11})

    for ax in axs[2, :]:
        ax.set_xlabel('Date')
        ax.set_xticklabels(pd.date_range(day_sdate, day_edate, freq='H'),
                           rotation=90)
        ax.xaxis.set_major_locator(dates.DayLocator())
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d %m'))
        ax.xaxis.set_minor_locator(dates.HourLocator(interval=3))
        ax.xaxis.set_minor_formatter(dates.DateFormatter('%H'))
        ax.tick_params(axis='x', which='major', pad=15)

    fig.tight_layout()
    fig.subplots_adjust(left=0.13, bottom=0.3)

    h1, l1 = axs[2, 2].get_legend_handles_labels()

    axs[2, 0].legend(h1, ('$\mathregular{SapEast_{sim}}$', '$\mathregular{SapWest_{sim}}$',
                          '$\mathregular{SapEast_{obs}}$', '$\mathregular{SapWest_{obs}}$'),
                     frameon=True, bbox_to_anchor=(-0., -1.5, 2, .102),
                     loc=3, ncol=8)

    fig.savefig('fig_14.png')

    pyplot.show(fig)


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

    fig, axs = pyplot.subplots(nrows=2, ncols=2, sharey='row', sharex='all')
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
        ax.text(0.05, 0.9, ('(a)', '(b)', '(c)', '(d)')[i], transform=ax.transAxes)

    axs[0, 0].set(
        ylabel='$\mathregular{A_{n, plant}\/[\mu mol\/s^{-1}]}$',
        ylim=(-20, 50))
    axs[1, 0].set(
        xlabel='Date',
        ylabel='$\mathregular{E_{plant}\/[g\/h^{-1}]}$',
        xlim=(beg_date, end_date),
        ylim=(-200, 1600))
    axs[1, 1].set(xlabel='Date')
    [ax.set_ylim(0, 5) for ax in fig.get_axes()[-2:]]
    fig.get_axes()[-1].set(ylabel='VPD [kPa]')

    for ax in axs[1, :].flatten():
        ax.set_xticklabels(datet, rotation=90)
        ax.xaxis.set_major_locator(dates.DayLocator())
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d %m'))

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.25)

    h1, l1 = axs[1, 1].get_legend_handles_labels()
    h2, l2 = fig.get_axes()[-1].get_legend_handles_labels()

    axs[1, 1].legend(h1 + h2, ['obs'] + ['sim%d' % d for d in range(5)] + ['VPD'], frameon=True,
                     bbox_to_anchor=(-1.5, -0.7, 2, .102), loc=3, ncol=8,
                     prop={'size': 11})

    fig.savefig('fig_15.png')

    pyplot.show(fig)


if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    from pathlib import Path
    from datetime import datetime, timedelta
    from scipy.stats import linregress
    from matplotlib import dates, pyplot

    from hydroshoot.architecture import mtg_load
    from hydroshoot.utilities import vapor_pressure_deficit as VPDa

    pyplot.style.use('seaborn-ticks')

    example_pth = Path(__file__).parents[2] / 'example'

    plot_figure_11()
    plot_figure_12()
    plot_figure_13()
    plot_figure_14()
    plot_figure_15()
    # # plot_meteo_vsp()
