#!/usr/bin/env python2
# -*- coding: utf-8 -*-


def plot_figure_13():
    """Generates figure 13 of the paper
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

        obs_grouped = obs_tab.groupby('itime').aggregate(mean)
        obs_grouped['itime'] = obs_grouped.index

        m_df = pd.merge(obs_grouped, sims)

        axs[0, 2].plot(m_df['An_plante'], m_df['An'], ('bo', 'ro')[iax])
        axs[1, 2].plot(m_df['E_plante'], m_df['E'], ('bo', 'ro')[iax])

        for egas, igas in enumerate(('An', 'E')):
            x = m_df[igas + '_plante'].values
            y = m_df[igas].values

            xindex = isfinite(x)
            yindex = isfinite(y)
            xyindex = xindex * yindex
            x, y = x[xyindex], y[xyindex]

            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            r2 = round(r_value ** 2, 3)
            bias = (y - x).mean()
            rmse = sqrt(((x - y) ** 2).mean())

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
    #pyplot.close(fig)


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


if __name__ == '__main__':
    import pandas as pd
    from pathlib import Path
    from datetime import datetime
    from numpy import mean, isfinite, sqrt
    from scipy.stats import linregress

    from matplotlib import dates, pyplot

    pyplot.style.use('seaborn-ticks')

    example_pth = Path(__file__).parents[2] / 'example'

    plot_figure_13()
    plot_meteo_vsp()
