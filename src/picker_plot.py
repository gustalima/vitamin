import numpy as np
import matplotlib.pyplot as plt

from gi.repository import Gtk
import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar
import matplotlib.pyplot as plt


class PointBrowser:
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """
    def __init__(self):
        self.lastind = 0

        self.text = ax.text(0.05, 0.95, 'selected: none',
                            transform=ax.transAxes, va='top')
        self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.4,
                                  color='yellow', visible=False)

    def onpress(self, event):
        if self.lastind is None: return
        if event.key not in ('n', 'p'): return
        if event.key=='n': inc = 1
        else:  inc = -1


        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
        self.update()

    def onpick(self, event):

       if event.artist!=line: return True

       N = len(event.ind)
       if not N: return True

       # the click locations
       x = event.mouseevent.xdata
       y = event.mouseevent.ydata


       distances = np.hypot(x, y)
       indmin = distances.argmin()
       dataind = event.ind[indmin]

       self.lastind = dataind
       self.update()

    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

#        ax2.cla()
#        ax2.plot(X[dataind])

#        ax2.text(0.05, 0.9, 'mu=%1.3f\nsigma=%1.3f'%(xs[dataind], ys[dataind]),
#                 transform=ax2.transAxes, va='top')
#        ax2.set_ylim(-0.5, 1.5)
        self.selected.set_visible(True)
        self.selected.set_data(xs[dataind], ys[dataind])

        self.text.set_text('selected: %d'%dataind)
       
        ano = plt.annotate(
                'teste', 
                xy = (xs[dataind], ys[dataind]), xytext = (-20, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'black', alpha = 0.3),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
                

        fig.canvas.draw()


if __name__ == '__main__':


        phi, psi = [12,3,45,5,7,8,45,345,12,3,123],[12,3,4,25,47,86,485,39,10,93,13]
       
        win = Gtk.Window()
        win.connect("delete-event", Gtk.main_quit )
        win.set_default_size(950,900)
        win.set_title("Ramachandran Plots")

        
        xs = phi
        ys = psi

        fig, (ax) = plt.subplots(1)
        ax.set(aspect='equal')
        ax.axis ([-180,180,-180,180])
        ax.imshow (plt.imread('../gui/regions.png')   , aspect='equal', extent = (-180,180,-180,180))
        ax.plot (phi, psi,  'o')
        ax.set_xlabel(r'$\phi$', fontsize=40)
        ax.set_ylabel(r'$\psi$', fontsize=40)
        ax.set_title ('General\n', fontsize=20)

        ax.set_xticks(np.arange(-180,180+1, 30))
        ax.set_yticks(np.arange(-180,180+1, 30))
        line, = ax.plot(xs, ys, '.', picker=5)  # 5 points tolerance

        browser = PointBrowser()        

        vbox = Gtk.VBox()
        win.add(vbox)
        canvas = FigureCanvas(fig)  # a Gtk.DrawingArea
        vbox.pack_start(canvas, True, True, 0)

        # Create toolbar
        toolbar = NavigationToolbar(canvas, win)
        vbox.pack_start(toolbar, False, False, 0)

        fig.canvas.mpl_connect('pick_event', browser.onpick)
        fig.canvas.mpl_connect('key_press_event', browser.onpress)




        win.show_all()
        Gtk.main()





