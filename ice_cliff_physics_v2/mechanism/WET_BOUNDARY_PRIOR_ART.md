# Wet-boundary prior-art audit,1November2026

## Buri2016: direct primary-method inspection

The [full publisher article](https://doi.org/10.1002/2016JF004039), §4.1.1, defines a pond-influence zone using geometric buffers and adds a constant0.033m/day subaqueous melt component. Its §6.3 identifies prescribed water melt and unrepresented variable water levels as limitations. Reburial is treated through slope/debris-view criteria, not a transient stream–debris–ice Robin boundary. This supports comparing a more explicit local boundary against that limited prescription, but does not make water-controlled cliff evolution a new discovery. It also does not prove no later model has such a boundary.

## What the public CliffEBM release actually implements

Pinned repository: [pburi/CliffEBM](https://github.com/pburi/CliffEBM/tree/ea5ca8bae9b3c67ec3cf2a8b56c724d87bb1909c).

`CliffEBM.R` Git blob ab54bc7cfad93dbfbc0ce80b542e9842993603a6 (76,640bytes) and `Functions.R` blob5f8f9331b6a0134f62131c1ef895716fdc20f59d were read. The first is1900lines and the second385. R lines1174–1194 calculate longwave radiation from surrounding debris. Lines1298–1304 sum radiative/sensible/latent atmospheric fluxes into melt. Literal searches for pond/subaqueous return no matches in the main file, but the stronger evidence is the actual flux expression and release framing: this is the gridded atmospheric ice-cliff EBM example, not a demonstrated reproduction of the full2016 dynamic pond/reburial implementation. It cannot be labelled the full Buri2016 comparator.

No full third-party source files are redistributed by this audit; immutable URLs and blob IDs are enough to retrieve the inspected source.

## Related heat/water physics

- [Reid and Brock2010](https://www.cambridge.org/core/journals/journal-of-glaciology/article/an-energybalance-model-for-debriscovered-glaciers-including-heat-conduction-through-the-debris-layer/21379FE9ED505F092BCA58276652B7A5): dynamic debris conduction/storage are established.
- [Giese et al.2020](https://tc.copernicus.org/articles/14/1555/2020/), §§3.3.1,5.5: the model resolves debris moisture/thermal properties but omits advective heat within debris. Reduced melt under wet debris cannot be attributed to the present external stream-film mechanism merely because its sign matches; evaporation and property changes are distinct causes.
- [Winter-Billington et al.2022](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2022.796877/full): SHAW-Glacier explicitly includes heat advection by rainwater and lateral flow and compares against debris temperature/ablation observations. Thus a general first coupled water/debris heat model claim is untenable.
- [Collier et al.2014](https://tc.copernicus.org/articles/8/1429/2014/): moisture reservoirs, phase changes and moisture-related cooling have precedent. Evaporative cooling is not signed sensible-heat export in a stream.
- [Petersen et al.2024](https://esurf.copernicus.org/articles/12/727/2024/): debris evacuation/stream undercutting/abandonment are established. The paper reports turbid streams, making unchanged absorbed debris radiation under submergence an assumption to audit rather than a field fact.

A primary Franz Josef2012 publisher search result discusses runoff intercepting and advecting incoming energy, but the full article returned403 in this audit; no equation-level claim is made from that result. Therefore even the broader water-cooling insight cannot responsibly be labelled previously unexamined. A specific stream-film export criterion applied to cliff toe rewetting remains a bounded candidate application, with historical priority unresolved.
