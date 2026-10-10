#!/usr/bin/env python3
"""Short direct-density validation figure before long integration."""
import argparse
import json
from pathlib import Path
import os
import tempfile
import sys
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'sat_string_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import prepare_output_directory
from sat_string.nrlmsis import OrbitGeometry


def japanese_style():
    fonts={x.name for x in font_manager.fontManager.ttflist}
    for name in ['Yu Gothic','Meiryo','Noto Sans CJK JP','IPAexGothic']:
        if name in fonts:plt.rcParams['font.family']=name;break
    else:raise RuntimeError('Japanese font required')
    plt.rcParams.update({'axes.unicode_minus':False,'axes.grid':True,'grid.alpha':.25,'font.size':10})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=prepare_output_directory(args.output);japanese_style()
    with np.load(args.run/'results.npz') as d:
        t=d['time_s'];rho=d['density_kg_m3'];a=d['a_m'];u=d['lambda_rad']
    spec=json.loads((args.run/'density_environment.json').read_text(encoding='utf-8'))['environment']
    g=OrbitGeometry(**{k:spec[k] for k in ['epoch_utc','inclination_deg','raan_deg','initial_phase_deg']});lon,lat,alt=g.coordinates(t[:,None],a,u)
    lst=(t[:,None]/3600+lon/15)%24
    fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True)
    for sat in [1,6,11,16,21]:axs[0].plot(t/3600,rho[:,sat-1]*1e12,label=f'衛星{sat}',lw=1)
    axs[0].set(ylabel='質量密度 [×10^-12 kg/m³]',title='衛星別密度：全LOW・人工外乱なし');axs[0].legend(ncol=5)
    axs[1].plot(t/3600,lat[:,10],label='中央衛星の測地緯度');axs[1].set(ylabel='測地緯度 [度]');axs[1].legend()
    axs[2].scatter(t/3600,lst[:,10],s=2,label='中央衛星の地方太陽時');axs[2].set(ylabel='地方太陽時 [時]',xlabel='開始からの時間 [時間]',ylim=(0,24));axs[2].legend()
    fig.suptitle('NRLMSIS 2.0の直接評価：2024年1月1日UTC開始・約3軌道\n地方太陽時は幾何学的な食判定ではありません')
    fig.tight_layout();fig.savefig(out/'density_validation.png',dpi=150,bbox_inches='tight');plt.close(fig)


if __name__=='__main__':main()
