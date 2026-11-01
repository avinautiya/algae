"""Observed-station input audit, never a simulated cliff result."""
from collections import defaultdict
from pathlib import Path
from statistics import mean
from audit_meteorology import load

def plot(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from datetime import date
    groups=defaultdict(list)
    for r in load(root/'meteodata_AWS2013.txt'):groups[r['timestamp_local_label'][:10]].append(r)
    # Do not present the last single midnight observation as a daily mean.
    days=sorted(d for d,rows in groups.items() if len(rows)==24)
    x=[date.fromisoformat(d) for d in days]
    def series(k):
        return [mean(r['values'][k] for r in groups[d] if r['values'][k] is not None) for d in days]
    fig,ax=plt.subplots(2,1,figsize=(10,6),sharex=True)
    ax[0].plot(x,series('TA'),color='#075985',label='Daily mean observed air temperature')
    ax[0].set_ylabel('Air temperature (°C)');ax[0].legend()
    ax[1].plot(x,series('SWIN'),label='Observed station shortwave')
    ax[1].plot(x,series('LWIN'),label='Observed station longwave')
    ax[1].set_ylabel('Daily mean radiation (W/m²)');ax[1].legend()
    ax[1].set_xlabel('Complete 24-hour date labels only; source timezone not established here')
    for a in ax:a.grid(alpha=.2)
    fig.suptitle('Observed Lirung AWS example from CliffEBM — input audit, not a cliff forecast')
    fig.tight_layout();fig.savefig(root/'observed_weather.png',dpi=180);fig.savefig(root/'observed_weather.pdf');plt.close(fig)
if __name__=='__main__':plot(Path(__file__).parent)
