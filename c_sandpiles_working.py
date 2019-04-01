'''
=====================================
Sand piles for comp - Will MCDONALD
======================================
'''

from random import randint, uniform
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sp

# for plotting
import matplotlib.pyplot as plt
from matplotlib import pylab
from numpy import linalg as LA
from scipy import stats



# --------------------------------------------------------------------------
class SandGrid:

    class Avalanche:

        def __init__(self):
            # Want AV data every time we have a AV
            self.time = 0
            self.av_at_time = []
            self.av_size = []
            self.av_lifet = []
            self.av_area = []
            self.av_r = []

        def nextdrop(self):
            self.time += 1

        def avalanche_occured(self, av_dict):
            # Store the step time of occurance
            self.av_at_time.append(self.time)

            # Store size of this averlance in list
            self.av_size.append(av_dict.get('size'))

            #Store area
            self.av_area.append(av_dict.get('area'))

            # Store lifetime of avalanche
            self.av_lifet.append(av_dict.get('lifetime'))

            # Store the radius of the averlance
            self.av_r.append(av_dict.get('radius'))




        def graph_distribution(self, array, title, xlabel, fit=False):
            # get your array, (should already be non-zero)
            vals = np.array(array)

            unique, counts = np.unique(vals, return_counts=True)

            if fit == True:
                x = np.log(unique)
                # Normalize the counts
                y = np.log(counts)#/np.array(np.linalg.norm(counts)))


                fig, ax = plt.subplots()

                # calculate polynomial
                z = np.polyfit(x, y,1)
                f = np.poly1d(z)
                print(f)

                # Print the relation on the graph
                plt.text(0.85,0.85, f,horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)

                # calculate new x's and y's
                x_new = np.linspace(min(x), max(x))
                y_new = f(x_new)
                plt.plot(x, y,x_new, y_new)
                ylabel = 'Frequency (log)'
                xlabel = xlabel + ' (log)'

                title = title + " (Normalized)"


            else:
                x = unique
                # Normalize the counts
                y = counts#/np.array(np.linalg.norm(counts))
                plt.plot(x, y)
                ylabel = 'Frequency'


            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.show()



        def plot_relationship(self, xlist, ylist, xlabel, ylabel, log=True):

            if log:
                x = np.log(np.array(xlist))
                y = np.log(np.array(ylist))
                scale = 'log'
            else:
                x = xlist
                y = ylist
                scale = 'lin'

            fig, ax = plt.subplots()

            # calculate polynomial
            z = np.polyfit(x, y,1)
            f = np.poly1d(z)
            print(f)

            # Print the relation on the graph
            plt.text(0.85,0.4, f,horizontalalignment='center', verticalalignment='center',transform=ax.transAxes)

            # calculate new x's and y's
            x_new = np.linspace(min(x), max(x))
            y_new = f(x_new)

            plt.plot(x, y,'o', x_new, y_new)

            title = 'Avalanche ' + xlabel + ' vs ' + ylabel + " (" + scale + scale +")"
            plt.title(title)
            plt.xlabel(xlabel + " (" + scale +")")
            plt.ylabel(ylabel + " (" + scale +")")

            ax = plt.gca()
            ax.set_facecolor((0.898, 0.898, 0.898))
            fig = plt.gcf()
            plt.show()



    # Initilze MxN grid
    # Outside is off the table
    def __init__(self, M, N, controls, add_sand=1, max=4, grad=45):
        self.M = M
        self.N = N
        self.controls = controls
        self.add_sand = add_sand
        self.max = max

        # Start with a 2D
        self.grid = [[0 for i in range(self.N)] for j in range(self.M)]
        self.mass_array = [0]
        # make a subclass of AV data
        self.avalanche = self.Avalanche()
        self.sand_diff = np.round(np.tan((grad * np.pi)/180)) - 1 # what the gradient corresponds to in sand difference.

        if self.sand_diff < 1:
            self.sand_diff = 1
        if self.sand_diff > 6:
            self.sand_diff = 6

        if self.controls.get('slope_table'):
            for i in range(self.M):
                for j in range(self.N):
                    self.grid[i][j] -= j

        self.slope_table = [[0 for i in range(self.N)] for j in range(self.M)]
        for i in range(self.M):
            for j in range(self.N):
                self.slope_table[i][j] -= j


    def run_sim(self):

        # Drop sand on Table
        # and return that location to .is_hill to see if this drop causes am avalanche
        # We return the location to calculate the radius
        av_dict = self.is_hill(self.drop_sand())
        # We have gone though an entire grid and there are no more spills
        # Check if we have had an avalanche at this time -> so we need to track if we did
        if  av_dict.get('size') != 0:

            # if we did, store the time it occured at
            self.avalanche.avalanche_occured(av_dict)

        # Now check the total number of sand on the table After everything has
        # spilled
        self.mass_array.append(self.grid_sum())



    # prints the sand value in a MxN grid
    def print_grid(self):

        for i in self.grid:
            for j in i:
                print("{:<2}".format(j),end=' ')
            print()
        print()

    # Randomly drop a peice of sand somewhere on Grid
    # Arrays start at 0
    def drop_sand(self):

        # every time we drop sand, increase time (AV) by 1
        self.avalanche.nextdrop()

        # if we are using control to only drop in middle
        if self.controls.get('drop_mid'):
            x = int(self.M/2)
            y = int(self.N/2)

        elif self.controls.get('drop_half'):
            x = randint(0, int(self.M/2))
            y = randint(0, self.N - 1)

        else:

            if self.controls.get('dis') == 'uniform':
                x = randint(0, self.M - 1)
                y = randint(0, self.N - 1)

            if self.controls.get('dis') == 'gauss':
                # want gauss centered on middle, with sigma such that being out off our plot happens 100 - 99.73% of the time (and distregard otherwise)
                # so self.M/2 to get the range either side, then /3 to be within 3 Stddiv
                x = -1
                y = -1

                while x not in range(0, self.M - 1):
                    x = int(np.random.normal(self.controls.get('meanx'), self.controls.get('sigmax')))

                while y not in range(0, self.N - 1):
                    y = int(np.random.normal(self.controls.get('meany'), self.controls.get('sigmay')))


            if self.controls.get('dis') == 'powerlaw':
                x = -1
                y = -1

                a = self.controls.get('powerlaw_a')
                while x not in range(0, self.M - 1):
                    x = int(np.random.power(a))

                while y not in range(0, self.N - 1):
                    y = int(np.random.power(a))

        self.grid[x][y] += self.add_sand

        return {'x': x, 'y':y}


    # Checks every location to see if it is too tall
    # can spill up hill...
    def is_hill(self, xy_dict):
        spill = True
        av = 0
        lifetime = 0

        # use this to track the locations of unique topping sites (area)
        tops_at = [[0 for i in range(self.N)] for j in range(self.M)]

        while spill:
            spill = False
            # check everywhere
            site_loc = []
            for i in range(self.M):
                for j in range(self.N):
                    # if is too big, spill

                    if self.grid[i][j] >= self.max:
                        site_loc.append([i,j])
                        # We have a site (Only count once)
                        tops_at[i][j] = 1

            for loc in site_loc:

                # avalanche adds the number of sand spilled during self.spill()
                if self.controls.get('smart_spill'):
                    spill_results = self.smart_spill(loc)
                    av += spill_results[0]
                    spill = spill_results[1]
                    # maybe we have a spill (might not have high enough gradient)
                else:
                    av += self.spill(loc)
                    #We have had a spill, so will need to check if there are any left over big piles.
                    spill = True

            for i in range(self.M):
                for j in range(self.N):
                    if self.grid[i][j] < self.slope_table[i][j]:
                        print('Gone below')
                        self.print_grid()
                        break
            if spill:
                # We have to increase the life time of avalanche by 1. But the previous
                lifetime += 1


        # Now we should count all the unique sites
        n_u_top = self.grid_sum(tops_at)

        # Claculate the radius from (x_0, y_0)
        max_r = self.radius_from(xy_dict, tops_at)

        # Return the size of avalanche in all
        return {'size': av , 'area':n_u_top, 'lifetime':lifetime, 'radius':max_r}


    def radius_from(self, xy_dict, site_grid):
        # We will consider the radius to be number of sites to the furtherest away site that is changed (Including spilling and stopping on)
        '''
        0  0  0  0
        0  0  1  1
        0  0  1  *
        0  0  0  1

        0  0  #  #
        0  #  1  1
        0  #  1  *
        0  0  #  1

        This is a 4x4 grid that started at * and 1 are unique topping sites. The radius of this is 2. (# is where sand spilled to but didnt cause a spill)
        '''
        # Site original drop
        x = xy_dict.get('x')
        y = xy_dict.get('y')
        top_loc = [[x,y]]

        # Get the locations in an array
        for i in range(self.M):
            for j in range(self.N):
                if site_grid[i][j] == 1:
                    top_loc.append([i,j])

        d_max = 0
        # Go though each location and find the linear linear distance
        sites_affected = []
        for a in top_loc:
            sites_affected.append(a)
            i = a[0]
            j = a[1]

            if i != 0:
                # Not on a hoz edge
                sites_affected.append([i-1,j])
            if i != self.M - 1:
                sites_affected.append([i+1,j])

            # Vert edge:
            if j != 0:
                sites_affected.append([i,j-1])
            if j != self.N - 1:
                sites_affected.append([i,j+1])

        #Now every point that was affected is in sites_affected
        for a in sites_affected:
            dx = abs(a[0] - x)
            dy = abs(a[1] - y)



            d_xy = dx + dy
            if d_max < d_xy:
                d_max = d_xy
        return d_max

    def spill(self, loc):
        i = loc[0]
        j = loc[1]
        # decrease the number of sand at (i,j) by self.max

        self.grid[i][j] -= self.max
        # check if on the edge
        # Can spill uphill - nonphysical

        if i != 0:
            # Not on a hoz edge
            self.grid[i-1][j] += 1
        if i != self.M - 1:
            self.grid[i+1][j] += 1

        # Vert edge:
        if j != 0:
            self.grid[i][j-1] += 1
        if j != self.N - 1:
            self.grid[i][j+1] += 1
        return self.max

    def smart_spill(self, loc):
        i = loc[0]
        j = loc[1]
        spill = False

        n_spilled = 0

        # there is a chance that the pile wont spill at all. (only spill stoc_prob percent of the time)
        if self.controls.get('stocastic_spill') and uniform(0, 1) > self.controls.get('stoc_prob'):
            return [n_spilled, spill]


        # check if on the edge
        # Can only spill downhill

        # Not on a hoz edge
        # Check to see if can spill downhill with a bigger gradient difference

        if i != 0:
            if self.grid[i][j]  > self.grid[i-1][j] + self.sand_diff:
                self.grid[i-1][j] += 1
                n_spilled += 1
        else:
            # spilled off edge
            n_spilled += 1

        if (i != self.M - 1):
            if self.grid[i][j] > self.grid[i+1][j] + self.sand_diff:
                self.grid[i+1][j] += 1
                n_spilled += 1
        else:
            # spilled off edge
            n_spilled += 1

        # Check if not on Vert edge and is spilling downhill:
        if j != 0:
            if self.grid[i][j] > self.grid[i][j-1] + self.sand_diff:
                self.grid[i][j-1] += 1
                n_spilled += 1
        else:
            n_spilled += 1

        if (j != self.N - 1):
            if self.grid[i][j] > self.grid[i][j+1] + self.sand_diff:
                self.grid[i][j+1] += 1
                n_spilled += 1
        else:
            n_spilled += 1


        self.grid[i][j] -= n_spilled


        # Note that is nothing spilled, it must be surrounded by hills over max, so running program again will get the surrounding first. It will evetually each somewhere it can spill/edge so will terminate.
        if n_spilled != 0:
            spill = True
        return [n_spilled, spill]

    def graph(self):

        # make 1d array for x and y direction ([0, 1, 2, 3...] etc)
        x = range(int(self.M))
        y = range(self.N)

        # get Z values from Table.grid - this is the num sand at a particular (x,y). Has to be a np array to work.
        Z = np.asarray(self.grid)

        fig, ax = plt.subplots()
        cax = ax.imshow(self.grid, cmap='terrain')
        ax.set_title('No slope - Random Placement, ({0},{1})'.format(self.M, self.N))
        cbar = fig.colorbar(cax)
        plt.show()
        #fig.savefig('./Graph.pdf')

    def grid_sum(self,grid=None):

        # map(sum, input) will return a list with the sums across all your rows,
        # then, the outer most sum will add up that list.
        if grid==None:
            return sum(map(sum, self.grid))
        else:
            return sum(map(sum, grid))

    def graph_array(self, Tables, xlabel, ylabel, title, steps):

        #plt.plot(arrays[1],'o',alpha=0.2, label="({0}x{1})".format(2*self.M, 2*self.N))
        plt.plot(arrays,'x', alpha=0.5, label="({0}x{1})".format(self.M, self.N) )

        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        title = '{0} ({1} steps)'.format(title, steps)
        plt.title(title)
        plt.legend(loc='upper right')
        plt.show()


def graph_array_size(Tables, xlabel, ylabel, title, steps):
    for Table in Tables:

        #plt.plot(arrays[1],'o',alpha=0.2, label="({0}x{1})".format(2*self.M, 2*self.N))
        plt.plot(Table.avalanche.av_size,'x', alpha=0.2, label="({0}x{1})".format(Table.M, Table.N) )

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    title = '{0} ({1} steps)'.format(title, steps)
    plt.title(title)
    plt.legend(loc='upper right')
    plt.show()


# --------------------------------------------------------------------------
if __name__ == '__main__':

    '''===================================================================='''
    # Controls
    M = 7
    N = 8
    steps = 20000
    gradient = 45 # 45 < grad < 80

    #Drop Distribution + turn on/off various parmas
    controls =  { 'dis'   : 'uniform',
                            #'gauss', # centered on middle
                            # 'powerlaw', # centered on left side
                 'meanx'  : M/2, # only used for gauss
                 'meany'  : N/2,
                 'sigmax' : M/6,
                 'sigmay' : N/6,
                 'powerlaw_a' : 1.1, # a > 1
                 'smart_spill' : False,
                 'stocastic_spill': False,
                 'stoc_prob' : 0.9,
                 'drop_mid' : False, # drop in the middle of the grid
                 'drop_half': False,# drop only on top half
                 'slope_table' : False
                 }


    '''===================================================================='''

    print('Running a {0}x{1} grid for {2} steps with a {3} Distribution \n'.format(M,N,steps, controls.get('dis')))
    # Initilze grid MxN
    # it will be flat z = 0
    Table = SandGrid(M, N, controls=controls, grad=gradient)
    Subtable = SandGrid(M, N, controls=controls, grad=gradient)



    # Run multiple ensembles:
    for ensembles in range(10):


        for i in range(steps):
            Subtable.run_sim()
            #Table2.run_sim()
    #Store data in list then input the average of all the runs into Table(master)



    '''
        # Drop sand on Table
        # and return that location to .is_hill to see if this drop causes am avalanche
        # We return the location to calculate the radius
        av_dict = Table.is_hill(Table.drop_sand())
        # We have gone though an entire grid and there are no more spills
        # Check if we have had an avalanche at this time -> so we need to track if we did
        if  av_dict.get('size') != 0:

            # if we did, store the time it occured at
            Table.avalanche.avalanche_occured(av_dict)

        # Now check the total number of sand on the table After everything has
        # spilled
        Table.mass_array.append(Table.grid_sum())'''
    '''
    Tables = [Table1, Table2]
    for Table in Tables:'''
    Table.print_grid()
    Table.graph()
    '''=============================================================='''
    #Table.graph_array([Table2.mass_array, Table1.mass_array], 'Time', 'Mass of Sand', 'Sand Mass Over Time',steps)

    '''Table.graph_array(Table.avalanche.av_size, 'Time', 'Avalanche Size', 'Avalanches Over Time',steps)
    Table.graph_array(Table.avalanche.av_area, 'Time', 'Number of Unique Steps', 'Avalanche Area over Time',steps)
    Table.graph_array(Table.avalanche.av_lifet, 'Time', 'Lifetime of Avalanche', 'Avalanche Lifetime over Time',steps)
    Table.graph_array(Table.avalanche.av_r, 'Time', 'Radius of Avalanche', 'Avalanche Radius over Time',steps)'''

    ''' ===================================================================='''
    #Table.avalanche.graph_distribution(Table.avalanche.av_size, 'Avalanche Size Distribution','Size',fit=False)
    Table.avalanche.graph_distribution(Table.avalanche.av_area, 'Avalanche Area Distribution', 'Area', fit=False)
    Table.avalanche.graph_distribution(Table.avalanche.av_lifet, 'Avalanche Lifetime Distribution', 'Lifetime', fit=False)
    Table.avalanche.graph_distribution(Table.avalanche.av_r, 'Avalanche Radius Distribution', 'Radius', fit=False)
    ''' ===================================================================='''

    Table.avalanche.plot_relationship(Table.avalanche.av_size, Table.avalanche.av_size, 'Size', 'Size')
    Table.avalanche.plot_relationship(Table.avalanche.av_area, Table.avalanche.av_size, 'Area', 'Size')
    Table.avalanche.plot_relationship(Table.avalanche.av_lifet, Table.avalanche.av_size, 'Life Time', 'Size')
    Table.avalanche.plot_relationship(Table.avalanche.av_r,Table.avalanche.av_size, 'Radius', 'Size')
