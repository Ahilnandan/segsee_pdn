"""Schematic of the modelled power delivery chain (pdn4 values), with the three droop loops marked."""
import schemdraw
import schemdraw.elements as e
import matplotlib; matplotlib.use('Agg')
schemdraw.config(fontsize=11, lw=1.6)
d = schemdraw.Drawing()
# --- off-chip buck
V = e.SourceV().up().label('3.3 V in', loc='left')
d += V
d += e.Line().right().length(0.5)
sw = e.Switch().right().label('buck switch\n2 MHz', loc='top')
d += sw
L1 = e.Inductor2().right().label('L1 2.2 µH\n40 mΩ', loc='top')
d += L1
nb = e.Dot()
d += nb
d.push(); d += e.Capacitor().down().label('board bulk\n47 µF\n(3 mΩ)', loc='bottom'); d += e.Ground(); d.pop()
# --- board path
d += e.Line().right().length(0.6)
d += e.Resistor().right().label('3 mΩ', loc='top')
d += e.Inductor2().right().label('board\n2 nH', loc='top')
npk = e.Dot()
d += npk
d.push(); d += e.Resistor().down().length(1.5).label('20 mΩ', loc='bottom'); d += e.Inductor2().down().length(1.5).label('100 pH', loc='bottom')
d += e.Capacitor().down().length(1.5).label('package\ndecap 1 µF', loc='bottom'); d += e.Ground(); d.pop()
# --- package / bumps
d += e.Line().right().length(0.6)
d += e.Resistor().right().label('50 mΩ', loc='top')
d += e.Inductor2().right().label('package/bumps\n0.5 nH', loc='top')
nd = e.Dot().label('DIE 1.1 V', loc='bottom', ofst=0.15)
d += nd
d.push(); d += e.Resistor().down().length(1.5).label('0.035 Ω', loc='bottom'); d += e.Capacitor().down().length(1.5).label('intrinsic\n3 nF', loc='bottom'); d += e.Ground(); d.pop()
d += e.Line().right().length(2.8)
d.push(); d += e.Resistor().down().length(1.5).label('damping\n0.2–0.7 Ω', loc='bottom'); d += e.Capacitor().down().length(1.5).label('added decap\n~3–11 nF', loc='bottom'); d += e.Ground(); d.pop()
d += e.Line().right().length(2.2)
d += e.SourceI().down().length(3).label('core load\ni(t) 0.1–0.5 A', loc='bottom'); d += e.Ground()
# loop annotations (plain text labels under the circuit)
fig = d.draw(show=False)
ax = fig.ax if hasattr(fig, 'ax') else fig.fig.axes[0]
y0 = min(ax.get_ylim()) - 0.6
for k, (txt, col) in enumerate([('1st droop: package/bumps 0.5 nH + die decap  ->  65-90 MHz, ~35 mV for a 150 mA step', '#eb6834'),
                                ('2nd droop: board 2 nH + package decap 1 uF  ->  ~3.6 MHz, ~5 mV', '#1baf7a'),
                                ('3rd droop: buck control loop + 47 uF board bulk  ->  ~100 kHz, ~5 mV', '#6250d6')]):
    ax.text(0, y0 - 0.9 * k, txt, color=col, fontsize=11, ha='left')
ax.set_ylim(y0 - 3.0, max(ax.get_ylim()))
fig.save('/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/pdn_circuit.png', dpi=140) if hasattr(fig, 'save') else d.save('/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/pdn_circuit.png', dpi=140)
print('saved')
