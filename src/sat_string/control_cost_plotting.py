"""Plots made from saved control-cost histories and tables only."""
from pathlib import Path
import numpy as np
from .plotting import plt
from matplotlib.colors import ListedColormap, BoundaryNorm


def plot_cost_report(cases, rows, config, output, dpi=200):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    names = ["all_low_drag_baseline", "controlled_disturbed", "controlled_undisturbed"]
    labels = ["A: disturbed, all LOW", "B: disturbed, controlled", "C: no disturbance, controlled"]
    n = config.formation.satellite_count
    paths = []

    def save(fig, name):
        fig.savefig(output/name, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths.append(name)

    cmap = ListedColormap(["#f1f5f9", "#f59e0b", "#172554", "#14b8a6"])
    norm = BoundaryNorm([-0.5,0.5,1.5,2.5,3.5], 4)
    fig, axes = plt.subplots(3, 1, figsize=(11,8), sharex=True, layout="constrained")
    for ax, name, label in zip(axes, names, labels):
        a = cases[name]
        h = a["time_s"]/3600
        im = ax.pcolormesh(h, np.arange(1,n+1), a["mode"].T, cmap=cmap, norm=norm, shading="nearest", rasterized=True)
        ax.set(title=label, ylabel="Satellite ID")
        ax.set_yticks([1,(n+1)//2,n])
    axes[-1].set_xlabel("Time (h)")
    bar = fig.colorbar(im, ax=list(axes), ticks=[0,1,2,3], shrink=0.85)
    bar.ax.set_yticklabels(["LOW","LOW to HIGH","HIGH","HIGH to LOW"])
    save(fig, "01_switching_states.png")

    for field, ylab, fname, offset in [
        ("gap_error_m","Spacing error (m)","02_spacing_histories.png",False),
        ("a_m","Change in a = change in altitude (m)","03_orbital_histories.png",True),
    ]:
        fig, axes = plt.subplots(3,1,figsize=(11,9),sharex=True,layout="constrained")
        colors = plt.colormaps["turbo"](np.linspace(0,1,cases[names[0]][field].shape[1]))
        for ax, name, label in zip(axes,names,labels):
            a = cases[name]
            values = a[field]-a[field][0] if offset else a[field]
            for i in range(values.shape[1]):
                ax.plot(a["time_s"]/3600,values[:,i],color=colors[i],lw=0.8)
            ax.set(title=label,ylabel=ylab)
            ax.grid(alpha=0.2)
            ax.ticklabel_format(axis="y",style="plain",useOffset=False)
            if not offset and name == "controlled_undisturbed":
                limit = max(float(np.max(np.abs(cases["controlled_disturbed"][field]))), 1e-6)
                ax.set_ylim(-1.1*limit, 1.1*limit)
                ax.text(0.02,0.88,f"max |e| = {np.max(np.abs(values)):.2e} m (roundoff)",
                        transform=ax.transAxes, fontsize=9)

        sm = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(1,len(colors)))
        fig.colorbar(sm,ax=list(axes),label="Satellite ID" if offset else "Gap ID",
                     ticks=np.unique(np.linspace(1,len(colors),5).astype(int)))
        axes[-1].set_xlabel("Time (h)")
        save(fig,fname)

    controlled = cases["controlled_disturbed"]
    low = cases["all_low_drag_baseline"]
    fig,ax = plt.subplots(figsize=(10,4.5),layout="constrained")
    ax.plot(controlled["time_s"]/3600, low["a_m"]-controlled["a_m"],lw=0.8)
    ax.set(xlabel="Time (h)",ylabel="Extra altitude loss B vs A (m)")
    ax.grid(alpha=0.2)
    save(fig,"03b_extra_orbital_loss.png")
    if not rows:
        return paths
    rows = sorted(rows,key=lambda r:r["h_on_m_s2"])
    x = np.array([r["h_on_m_s2"] for r in rows])/1e-7
    specs = [
        ("e_max_m","Maximum spacing error (m)","04_threshold_accuracy.png"),
        ("high_stable_total_h","Stable HIGH exposure (satellite h)","05_threshold_exposure.png"),
        ("extra_mean_a_loss_vs_low_m","Extra mean altitude loss vs A (m)","06_threshold_loss.png"),
        ("range_km","Maximum observed switching range (km)","07_threshold_range.png"),
        ("switch_count","Total transition starts","08_threshold_switch_count.png"),
    ]
    for key,ylab,name in specs:
        fig,ax=plt.subplots(figsize=(7,4.5),layout="constrained")
        ax.plot(x,[r[key] for r in rows],"o-",label="Controlled + disturbed")
        if key=="high_stable_total_h":
            ax.plot(x,[r["high_equivalent_total_h"] for r in rows],"s--",label="Area-equivalent HIGH (including slews)")
            ax.legend(fontsize=8)
        if key=="e_max_m":
            ax.axhline(np.max(np.abs(low["gap_error_m"])),color="gray",ls="--",label="A: all LOW")
            ax.legend()
        ax.axvline(config.controller.threshold_on_m_s2/1e-7,color="gray",alpha=0.45,ls=":")
        ax.set(xlabel="h_on (10^-7 m/s²); h_off = h_on - 10^-8 m/s²",ylabel=ylab)
        ax.grid(alpha=0.2)
        save(fig,name)

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    for ax,cost,ylab in [(axes[0],"extra_mean_a_loss_vs_low_m","Extra mean altitude loss (m)"),
                         (axes[1],"high_equivalent_total_h","Equivalent HIGH exposure (satellite h)")]:
        for j,r in enumerate(rows):
            ax.scatter(r["e_max_m"],r[cost],marker=["o","s","^","D","v","P","*"][j],
                       color=plt.colormaps["viridis"]((x[j]-x.min())/(x.max()-x.min())),
                       s=55,label=f"{r['h_on_m_s2']/1e-7:.3g}")
        ax.scatter([np.max(np.abs(low["gap_error_m"]))],[0],marker="x",color="black",label="A: all LOW")
        ax.set(xlabel="Maximum spacing error (m)",ylabel=ylab)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle("Accuracy–cost points; legend gives h_on / 10^-7 (lower left is preferable)")
    save(fig,"09_accuracy_cost.png")

    fig,ax=plt.subplots(figsize=(7,4.5),layout="constrained")
    for side,marker in [("left","o"),("right","s")]:
        y=[r[f"speed_{side}_m_s"] if r[f"speed_{side}_m_s"] is not None else np.nan for r in rows]
        ax.plot(x,y,marker+"--",label=side)
    ax.set(xlabel="h_on (10^-7 m/s²)",ylabel="Accepted first-front speed (m/s)",
           xlim=(x.min()-0.03,x.max()+0.03))
    ax.text(0.98,0.80,"Missing speeds: see fit status in case JSON",
            transform=ax.transAxes,ha="right",fontsize=8)
    ax.legend()
    ax.grid(alpha=0.2)
    save(fig,"10_threshold_front_speed.png")
    return paths



def plot_horizon_cost(cases, config, output, dpi=200):
    """Show cumulative area-equivalent exposure and orbital loss at 24/72 hours."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    low = cases["A_low"]
    labels = {"A_low":"A: all LOW", "B_current":"Current h_on=1.5e-7",
              "h_on_1.800e-07":"h_on=1.8e-7", "h_on_2.000e-07":"h_on=2.0e-7"}
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,layout="constrained")
    for name,arrays in cases.items():
        t=arrays["time_s"]
        weights=np.sum((arrays["area_m2"]-config.spacecraft.area_low_m2)/
                       (config.spacecraft.area_high_m2-config.spacecraft.area_low_m2),axis=1)
        cumulative=np.concatenate(([0.0],np.cumsum(np.diff(t)*(weights[1:]+weights[:-1])/2)))
        axes[0].plot(t/3600,cumulative/3600,label=labels[name],lw=1.6)
        axes[1].plot(t/3600,np.mean(low["a_m"]-arrays["a_m"],axis=1),label=labels[name],lw=1.6)
    axes[0].set_ylabel("Cumulative equivalent HIGH (satellite h)")
    axes[1].set(xlabel="Time (h)",ylabel="Extra mean altitude loss vs A (m)")
    for ax in axes:
        ax.axvline(config.simulation.duration_s/3600,color="gray",ls=":",label="Original 24 h endpoint")
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8)
    path=output/"11_horizon_cost.png"
    fig.savefig(path,dpi=dpi,bbox_inches="tight")
    plt.close(fig)
    return path
