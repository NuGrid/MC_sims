import utils
from matplotlib import pylab as plt
import numpy as np
import nugridse
import ppn
import astronomy as ast
from scipy import integrate
import operator

def plot_multi_species(p,species,limits,plot_title):
    i = 0
    for this_species in species:
        i += 1; a,b = utils.linestyle(i)
        p.plot('mass',this_species,logy=True,shape=a,markevery=b,limits=limits)
    title(plot_title)
    legend(loc=0)

def select_fig(ifig,close_fig=True):
    if close_fig:
        plt.close(ifig)
    plt.figure(ifig)

#CR: The following function allows to plot the simulation data from
#the Sakurai11 paper.
def plot_RUN48_sakurai_abund(datfile,zrange=None,plot_elements=[],msize=12,marker='o',color='b',linestyle='-',sol_ab='',decay=True,gridon=True,label=''):
    '''
    read in Asplund's 1999 Sakurai's element abundances (12.0 for H)

    Parameters 
    ----------
    datfile : string
       directory path to where the file of type ABUPP00779913000.DAT resides

    zrange : tuple, float, optional 
       range in Z for plot

   grid : boolean
        if true, plot grid  as in herwig11


    decay : boolean
       if True, decays abundance given by datfile

    sol_ab : string
       define the path to the solar abundance file you want to use for normalization, 
        e.g. iniab1.0E-02.ppn_asplund05

    '''

    #CR added to read RUN48
    import read_ABUPP as abupp
    reload(abupp)
    import utils as u
    import nugridse as mp
    from utils import get_z_from_el

    mass,t9ppg,rhoppg,dppg,rppg,dppgLagr,isotopes_massfrac,isotopes = abupp.read_abupp(file1=datfile)
    #in RUN48 case the following isotopes are available: NEUT     9 PROT    10 H   2   11 HE  3   12 HE  4   13 BE  7   14 B   8   15 LI  7   16 C  11   17 B  11   18 C  12   19 C  13   20 N  13

    #get all available elements from .DAT file
    elements=[]
    Zs=[]
    massfrac=[]
    for k in range(len(isotopes)):
          ele=isotopes[k].split('-')[0]
          if k==0:
                Zs.append(0)
                elements.append('Neutron') #as in SE files
                continue
          if ele not in elements:
                elements.append(ele)
                Zs.append(get_z_from_el(ele.capitalize()))

    # read in solar abundances in the path sol_ab
    if True:
        utils.solar(sol_ab,1.)
        sol_abu1=utils.solar_abundance
        sol_abu=[]
        names_sol=utils.names_sol
        for k in range(len(sol_abu1)):
                sol_abu.append(sol_abu1[names_sol[k]])
        n_sol=len(sol_abu)
        z_sol=utils.z_sol
        sol_abu_ele=[]
        z_sol_ele=[]
        name_sol_ele=[]
        import re
        for k in range(n_sol):
                name=re.split('(\d+)',names_sol[k])[0].capitalize().strip()
                if not z_sol[k] in z_sol_ele:
                        z_sol_ele.append(z_sol[k])
                        sol_abu_ele.append(sol_abu[k])
                        name_sol_ele.append(name)
                else:
                        idx=z_sol_ele.index(z_sol[k])
                        sol_abu_ele[idx] = sol_abu_ele[idx] + sol_abu[k]
    if False:
        import utils as u
        iniabu=u.iniabu(sol_ab)
        iniiso=[]
        iniabu_massfrac=[]
        for k in range(len(iniabu.habu)):
                iso=iniabu.habu.keys()[k]
                name=re.split('(\d+)',iso)[0].strip().capitalize()
                iniabu_massfrac.append(iniabu.habu.values()[k])
    yields=[0]*len(elements)

    #get surface abundances from .DAT files
    isotopes_massfrac_surf=[]
    for h in range(len(isotopes_massfrac)):
        isotopes_massfrac_surf.append(isotopes_massfrac[h][-1])

    #How to decay the abundance?
    #similar toaverage_iso_abund_marco, need for that the surf.h5 file

    import nugridse as mp
    sefiles=mp.se('/home/dpa/iprocess-multizone','M3.00Z.0001.0047801.surf.h5')
    isotopes1=sefiles.se.isotopes

    # assign isotopes to mass_frac vector, set rest to 0 
    # should be 2-D array, but populate only first entry
    sefiles.mass_frac=[[0]*len(isotopes1)]
    for k in range(len(isotopes1)):
        if isotopes1[k] in isotopes:
                idx=isotopes.index(isotopes1[k])
                sefiles.mass_frac[0][k]=isotopes_massfrac_surf[idx]

    #needs spe array, part of read_iso_abund_marco
    if decay:
        isotope_names = sefiles.se.isotopes
        u.convert_specie_naming_from_h5_to_ppn(isotope_names)
        names_ppn_world = u.spe
        number_names_ppn_world = u.n_array
        u.define_zip_index_for_species(names_ppn_world,number_names_ppn_world)
        u.stable_specie()
        stable_isotopes=u.stable #available stable isotopes due to decay
        sefiles.decay(sefiles.mass_frac)
        #from decay()
        #get decayed abundance, as chosen above get tonly first entry
        mass_frac_decayed_stable=mp.decayed_multi_d[0]
        mass_frac_decayed=[] #contains all isotopes, including unstable
        #print 'show decayed isotopes'
        for k in range(len(isotopes)):
                other_name_scheme=isotopes[k].split("-")[0].upper()+(5-len(isotopes[k])+1)*" "+isotopes[k].split("-")[1]
                if other_name_scheme in stable_isotopes:
                        idx=stable_isotopes.index(other_name_scheme)
                        mass_frac_decayed.append(mass_frac_decayed_stable[idx])
                else:
                        #unstable
                        mass_frac_decayed.append(1e-99)

    #for k in range(len(mass_frac_decayed)):
        #print isotopes[k],mass_frac_decayed[k]

    else: #no decay
        mass_frac_decayed=isotopes_massfrac_surf

    # convert isotopes in mass_frac_decayed into elements
    yields=[0]*len(elements)
    #print 'isotopes massfrac ',isotopes_massfrac_surf
    for iso in isotopes:
         ele=iso.split('-')[0]
         if ele in elements:
               yields[elements.index(ele)]= np.array(yields[elements.index(ele)]) + np.array( mass_frac_decayed[isotopes.index(iso)])


    #### Convert into spectroscopic notation
    #yields=yields[-1]
   # print 'take yields only at surface',mass[-1] 
    #print 'yields : ',yields
    yields_spec=[]
    Zs_spec=[]
    for k in range(len(elements)):
        if not elements[k] in plot_elements:
                continue
        if elements[k] in name_sol_ele:
                idx=name_sol_ele.index(elements[k])
                #if abundance too low skip
                if (yields[k])<1e-30:
                        print ('Warning: Isotope ',elements[k],yields[k],', skip element')
                        continue
                #print Zs[k],elements[k],yields[k],sol_abu_ele[idx]
                yields_spec.append(np.log10(yields[k]/sol_abu_ele[idx]))
                Zs_spec.append(Zs[k])
        else:
                print ('element ',elements[k],' not found., skip')
    #print 'final Zs_spec',Zs_spec
    #print yields_spec
    plt.plot(Zs_spec,yields_spec,marker=marker,color=color,linestyle=linestyle,markersize=msize,label=label)
    if zrange is not None:
        plt.xlim(zrange[0],zrange[1])
    plt.ylabel('[$X_i/X_{i,sol}$]')
    plt.xlabel('charge number Z')
    plt.minorticks_on()
    if gridon==True:
        plt.grid(b=True,which='major',color='k',linestyle=':')

# The following provides plotting of the Asplund Sakurai observation.
# This is an example of a repeatedly occuring task, i.e. access and
# plot observational abundance distribution data. There are two options:
# 1. use the Laurent observational database tool (FH would prefer this is
#    explored)
# 2. we use ascii_table module with the existing (or if needed new) table
#    write-read methods (less prefered)

# for now, let's just turn this functionality the way it is into a
# a method we can reuse in a notebook

def plot_johnson2004_abund(data_dir,zrange=None,yrange=None,symbol='bo',msize=12):
    '''
    read in Johnson's 2004 CS 31062-50 element abundances [X/Fe]

    Parameters 
    ----------
    data_dir : string
       directory path to where the file Johnson2004.txt is

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    symbol : string
       symbol to plot observational data

    msizel : integer
       size of the plotting symbol      
    '''
	
    lines=open(data_dir+'Johnson2004.txt','r').readlines()
    nelsak=np.size(lines)
    zsak=np.linspace(0,0,nelsak)
    elserryk=np.linspace(0,0,nelsak)
    elserrxk=np.linspace(0,0,nelsak)
    elname=["  " for x in range(nelsak)]
    elsak=np.linspace(0,0,nelsak)
    i=0
    for line in lines:
	parts=line.split()
	elname[i]=parts[0]
 	zsak[i]=np.array(float(parts[1]))
 	elsak[i]=np.array(float(parts[2]))
 	elserryk[i]=np.array(float(parts[3]))
	i += 1

    plt.errorbar(zsak,elsak,elserryk,elserrxk,symbol,mfc=None,markersize=msize)
    if zrange is not None:
        plt.xlim(zrange[0],zrange[1])
    if yrange is not None:
        plt.ylim(yrange[0],yrange[1])

def plot_asplund99_sakurai_abund(data_dir,zrange=None,yrange=None,msize=12,dontplotMay=True):
    '''
    read in Asplund's 1999 Sakurai's element abundances (12.0 for H)

    Parameters 
    ----------
    data_dir : string
       directory path to where the files Asplund1999_Sakurai.txt 
       and Asplund1999_Sakurai_May.txt are 

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    dontplotMay : boolean, optional
       when it is False, plot also observations taken in May
    '''
	
    lines=open(data_dir+'Asplund1999_Sakurai.txt','r').readlines()
    nelsak=np.size(lines)
    zsak=np.linspace(0,0,nelsak)
    elsol=np.linspace(0,0,nelsak)
    elname=["  " for x in range(nelsak)]
    elsak=np.linspace(0,0,nelsak)
    i=0
    for line in lines:
	parts=line.split()
 	zsak[i]=np.array(float(parts[0]))
 	elsol[i]=np.array(float(parts[1]))
	elname[i]=parts[2]
 	elsak[i]=np.array(float(parts[3]))
	i += 1
    for i in range(nelsak):
        if zsak[i] == 26:
            isakfe = i
    y=np.linspace(0,0,nelsak)
    feh=elsak[isakfe]-elsol[isakfe]
    for i in range(nelsak):
#	y[i]=elsak[i]-elsol[i]-feh
	y[i]=elsak[i]-elsol[i]
    xerr=np.linspace(0,0,nelsak)
    yerr=np.linspace(0,0,nelsak)
    for i in range(nelsak):
   	yerr[i]=0.3 
#   plt.plot(zsak,y,'k*',markersize=16)
    plt.errorbar(zsak,y,yerr,xerr,'ko',markersize=msize)
    #plt.plot(zsak,y,'k--')
    if zrange is not None:
        plt.xlim(zrange[0],zrange[1])
    if yrange is not None:
        plt.ylim(yrange[0],yrange[1])

    if dontplotMay:
	return nelsak, zsak, y, yerr # don't plot Sakuari's object May data 
    
    lines=open(data_dir+'Asplund1999_Sakurai_May.txt','r').readlines()
    nelsak=np.size(lines)
    zsak=np.linspace(0,0,nelsak)
    elsol=np.linspace(0,0,nelsak)
    elname=["  " for x in range(nelsak)]
    elsak=np.linspace(0,0,nelsak)
    i=0
    for line in lines:
	parts=line.split()
 	zsak[i]=np.array(float(parts[0]))
 	elsol[i]=np.array(float(parts[1]))
	elname[i]=parts[2]
 	elsak[i]=np.array(float(parts[3]))
	i += 1
    for i in range(nelsak):
        if zsak[i] == 26:
            isakfe = i
    y=np.linspace(0,0,nelsak)
    feh=elsak[isakfe]-elsol[isakfe]
    for i in range(nelsak):
#	y[i]=elsak[i]-elsol[i]-feh
	y[i]=elsak[i]-elsol[i]
    xerr=np.linspace(0,0,nelsak)
    yerr=np.linspace(0,0,nelsak)
    for i in range(nelsak):
   	yerr[i]=0.3 
#   plt.plot(zsak,y,'k*',markersize=16)
    plt.errorbar(zsak,y,yerr,xerr,'ko',mfc=None,markersize=msize)
    #plt.plot(zsak,y,'k--')
    if zrange is not None:
        plt.xlim(zrange[0],zrange[1])
    if yrange is not None:
        plt.ylim(yrange[0],yrange[1])

# The following provides plotting of first derivatives of the Asplund Sakurai observation.
# This is an example of a repeatedly occuring task, i.e. access and
# plot observational abundance distribution data. There are two options:
# 1. use the Laurent observational database tool (FH would prefer this is
#    explored)
# 2. we use ascii_table module with the existing (or if needed new) table
#    write-read methods (less prefered)

# for now, let's just turn this functionality the way it is into a
# a method we can reuse in a notebook
def plot_deriv_asplund99_sakurai_abund(data_dir,zrange=None,msize=12):
    '''
    read in Asplund's 1999 Sakurai's element abundances (12.0 for H)

    Parameters 
    ----------
    data_dir : string
       directory path to where the file Asplund1999_Sakurai.txt is

    zrange : tuple, float, optional 
       range in Z for plot
    '''
    
    lines=open(data_dir+'Asplund1999_Sakurai.txt','r').readlines()
    nelsak=np.size(lines)
    zsak=np.linspace(0,0,nelsak)
    elsol=np.linspace(0,0,nelsak)
    elname=["  " for x in range(nelsak)]
    elsak=np.linspace(0,0,nelsak)
    i=0
    for line in lines:
	parts=line.split()
 	zsak[i]=np.array(float(parts[0]))
 	elsol[i]=np.array(float(parts[1]))
	elname[i]=parts[2]
 	elsak[i]=np.array(float(parts[3]))
	i += 1
    for i in range(nelsak):
        if zsak[i] == 26:
            isakfe = i
    y=np.linspace(0,0,nelsak)
    feh=elsak[isakfe]-elsol[isakfe]
    for i in range(nelsak):
#	y[i]=elsak[i]-elsol[i]-feh
	y[i]=elsak[i]-elsol[i]
    for i in range(nelsak):
	if zsak[i] == 37 or zsak[i] == 38 or zsak[i] == 39:
		y[i]=y[i+1]-y[i]
	else:
		y[i]=-99
    xerr=np.linspace(0,0,nelsak)
    yerr=np.linspace(0,0,nelsak)
    for i in range(nelsak):
   	yerr[i]=0.3 
#   plt.plot(zsak,y,'k*',markersize=16)
    plt.errorbar(zsak,y,yerr,xerr,'ko',markersize=msize)
    if zrange is not None:
        plt.xlim(zrange[0],zrange[1])

# The following function compared decayed surface element  
# abundances for a same specified model from two specified mppnp H5_surf files
def comp_surf_abund(dir,comp_case,case,model,zrange=None,yrange=None,max_ratio=99.,decay=True,shape='',msize=14,labels=True,label_fsize=12,dely0=0.2,extr0=0.0,title=""):
    '''
    Parameters 
    ----------
    dir : string
       directory path to where H5_surf files are

    comp_case : string
       H5_surf_case is the directory in the dir with surf files to compare with

    case : string
       H5_surf_case is the directory in the dir with surf files to plot

    model : integer
       model number to compare

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    max_ratio : float
       if max(y_work) > max_ratio then print comp_case

    decay : boolean
       if True then decayed surface abundance are plotted

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any              
    '''

# read in model surface abundances
    s=nugridse.se(dir+"/H5_surf_"+case)
    scomp=nugridse.se(dir+"/H5_surf_"+comp_case)

    if decay:
    	el_abu_hif=s.get(model,'elem_massf_decay')
    	el_abu_ref=scomp.get(model,'elem_massf_decay')
    else:
    	el_abu_hif=s.get(model,'elem_massf')
    	el_abu_ref=scomp.get(model,'elem_massf')

    n_el_hif=len(el_abu_hif)
    z_hif=np.linspace(0,0,n_el_hif)
    y_hif=np.linspace(0,0,n_el_hif)
    el_name_hif=["  " for x in range(n_el_hif)]
    z_work=np.linspace(0,0,n_el_hif)
    y_work=np.linspace(0,0,n_el_hif)
    el_name_work=["  " for x in range(n_el_hif)]

    el_name_hif[0]='n'
    for i in range(n_el_hif):
    	z_hif[i]=float(i)       #  Z=i in mppnp surf data output
    	if (i>0):
        	el_name_hif[i]=utils.get_el_from_z(i)

# prepare plot arrays
    j=-1
    for i in range(n_el_hif):
	z_el=z_hif[i]
	if z_el != 0 and z_el != 43 and z_el != 61:
		j=j+1
        	z_work[j]=z_hif[i]
            	y_work[j]=-20.
		el_name_work[j]=el_name_hif[i]
            	if el_abu_ref[i]>1e-30 and el_abu_hif[i]>1e-30:
               		y_work[j]=np.log10(el_abu_hif[i]/el_abu_ref[i])
    n_plot=j

    z_plot=np.linspace(0,0,n_plot)
    y_plot=np.linspace(0,0,n_plot)
    el_name_plot=["  " for x in range(n_plot)]
    for i in range(n_plot):
    	z_plot[i]=z_work[i]
    	y_plot[i]=y_work[i]
	el_name_plot[i]=el_name_work[i]

    file_out='max_ratios.txt'
    fout=open(file_out,'a')

    max_index=zrange[0]+np.argmax(abs(y_work[zrange[0]:zrange[-1]]))
    max_y=max(abs(y_work[i]) for i in range(zrange[0],zrange[-1]))
    if max_y > max_ratio:
        fout.write(el_name_work[max_index]+"  "+str(max_y)+"  "+case+"\n")

# plot ratios of surface mppnp hif abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':   
    	shape='k-'     
    plt.plot(z_plot,y_plot,shape,markersize=msize)
    #plt.plot(z_plot,y_plot,'k-')
    if zrange == None:
    	plt.xlim(0,85)
    else:
    	plt.xlim(zrange)
    if yrange == None:
    	plt.ylim(-5,5)
    else:
    	plt.ylim(yrange)
    #plt.xlabel('Z')
    plt.xlabel('charge number Z')
    plt.ylabel('$\log_{10}\,X/X_\mathrm{RUN48}$')
    if labels == True:
    	for i in range(n_plot):
		if zrange != None and z_plot[i] >= zrange[0] and z_plot[i] <= zrange[-1]:
			dely = dely0
			extr = 0.0
    			if i%2 != 0:
				extr = extr0
			dely = dely + extr
			if i>0 and i<n_plot-2 and y_plot[i] < y_plot[i-1] and y_plot[i] < y_plot[i+1]:
        			dely = -0.2 - extr
			plt.text(z_plot[i],y_plot[i]+dely,el_name_plot[i],horizontalalignment='center',fontsize=label_fsize)
    plt.title(title)

# The following function plots solar-scaled decayed surface element  
# abundances for a specified model from a specified mppnp H5_surf files
def plot_surf_abund(sol_ab,dir,case,model,ref=-1,zrange=None,yrange=None,decay=True,clr='',mrk='o',mrkclr='k',shape='',msize=14,labels=True,label_fsize=12,dely0=0.2,extr0=0.0,title=""):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    ref : integer
       if ref > 0 then ref = model number to scale plotted abundances
       if ref = -1 then use sol_ab for the scaling
       if ref = -2 then still use sol_ab which are now in fact initial ones for the scaling 
       if ref = -10 then do not use any scaling

    dir : string
       directory path to where H5_surf files are

    case : string
       H5_surf_case is the directory in the dir with surf files to plot

    model : integer
       model number to plot

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    decay : boolean
       if True then decayed surface abundance are plotted

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any              
    '''

# read in solar abundances in the path sol_ab
    if ref == -1 or ref == -2:
    	utils.solar(sol_ab,1.)  
    	sol_abu=utils.solar_elem_abund
    	n_sol=len(sol_abu)

# read in model surface abundances
    if case == "":
    	s=nugridse.se(dir+"/H5_surf")
    #	sout=nugridse.se(dir+"/H5_out")
    else:
	s=nugridse.se(dir+"/H5_surf_"+case)
#	sout=nugridse.se(dir+"/H5_out_"+case)

    if decay:
    	el_abu_hif=s.get(model,'elem_massf_decay')
    else:
    	el_abu_hif=s.get(model,'elem_massf')

#   mass_hif=sout.get(model,'mass')
#   xli7_hif=sout.get(model,'Li-7')
#   xbe7_hif=sout.get(model,'Be-7')
#   shell_mass=0.
#   mass_li7=0.
#   mass_be7=0.
#   for ili7 in range(1,len(mass_hif)):
#       mass_li7=mass_li7+0.5*(xli7_hif[ili7]+xli7_hif[ili7-1])*(mass_hif[ili7]-mass_hif[ili7-1])
#       mass_be7=mass_be7+0.5*(xbe7_hif[ili7]+xbe7_hif[ili7-1])*(mass_hif[ili7]-mass_hif[ili7-1])
#   shell_mass=mass_hif[-1]-mass_hif[0]
#   xli7_shell=mass_li7/shell_mass
#   xbe7_shell=mass_be7/shell_mass
#   print "\n xli7_shell=",xli7_shell, "\n xbe7_shell=",xbe7_shell

    n_el_hif=len(el_abu_hif)
    z_hif=np.linspace(0,0,n_el_hif)
    y_hif=np.linspace(0,0,n_el_hif)
    el_name_hif=["  " for x in range(n_el_hif)]
    z_work=np.linspace(0,0,n_el_hif)
    y_work=np.linspace(0,0,n_el_hif)
    el_name_work=["  " for x in range(n_el_hif)]

    el_name_hif[0]='n'
    for i in range(n_el_hif):
    	z_hif[i]=float(i)       #  Z=i in mppnp surf data output
    	if (i>0):
        	el_name_hif[i]=utils.get_el_from_z(i)

    if ref >= 0:
	if decay:
    		el_abu_ref=s.get(ref,'elem_massf_decay')
	else:
    		el_abu_ref=s.get(ref,'elem_massf')
    	n_el_ref=len(el_abu_ref)
    	z_ref=np.linspace(0,0,n_el_ref)
    	y_ref=np.linspace(0,0,n_el_ref)
    	el_name_ref=["  " for x in range(n_el_ref)]

    	el_name_ref[0]='n'
    	for i in range(n_el_ref):
    		z_ref[i]=float(i)       #  Z=i in mppnp surf data output
    		if (i>0):
        		el_name_ref[i]=utils.get_el_from_z(i)

# prepare plot arrays
    if ref == -1 or ref == -2:
    	j=-1
    	for i in range(n_el_hif):
    		for k in range(n_sol):
        		z_sol=k+1
        		if float(z_sol)==z_hif[i] and z_sol != 43 and z_sol != 61:
            			j=j+1
            			z_work[j]=z_hif[i]
            			y_work[j]=-20.
				el_name_work[j]=el_name_hif[i]
            			if sol_abu[k]>1e-30 and el_abu_hif[i]>1e-30:
                			y_work[j]=np.log10(el_abu_hif[i]/sol_abu[k])
	#			if z_sol==3:
	#				print (xli7_shell+xbe7_shell),sol_abu[k]
	#				y_work[j]=np.log10((xli7_shell+xbe7_shell)/sol_abu[k])
    	n_plot=j+1
    if ref >= 0:
    	j=-1
    	for i in range(n_el_hif):
		z_el=z_hif[i]
		if z_el != 0 and z_el != 43 and z_el != 61:
			j=j+1
        		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if el_abu_ref[i]>1e-30 and el_abu_hif[i]>1e-30:
               			y_work[j]=np.log10(el_abu_hif[i]/el_abu_ref[i])
	#		if z_el==3:
	#			y_work[j]=np.log10((xli7_shell+xbe7_shell)/sol_abu[k])
	n_plot=j
    if ref == -10 :
    	j=-1
    	for i in range(n_el_hif):
		z_el=z_hif[i]
		if z_el != 0 and z_el != 43 and z_el != 61:
			j=j+1
        		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if el_abu_hif[i]>1e-30:
               			y_work[j]=np.log10(el_abu_hif[i])
	#		if z_el==3:
	#			y_work[j]=np.log10((xli7_shell+xbe7_shell)/sol_abu[k])
	n_plot=j

    z_plot=np.linspace(0,0,n_plot)
    y_plot=np.linspace(0,0,n_plot)
    el_name_plot=["  " for x in range(n_plot)]
    for i in range(n_plot):
    	z_plot[i]=z_work[i]
    	y_plot[i]=y_work[i]
	el_name_plot[i]=el_name_work[i]

# plot solar scaled surface mppnp hif abundances
    #font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    #plt.rc('font', **font)

    if shape == '':   
    	plt.plot(z_plot,y_plot,color=clr,marker=mrk,markerfacecolor=mrkclr,markersize=msize)
    else:
    	plt.plot(z_plot,y_plot,shape,markersize=msize)
    #plt.plot(z_plot,y_plot,'k-')
    if zrange == None:
    	plt.xlim(0,85)
    else:
    	plt.xlim(zrange)
    if yrange == None:
    	plt.ylim(-5,5)
    else:
    	plt.ylim(yrange)
    #plt.xlabel('Z')
    plt.xlabel('charge number Z')
    if ref == -1:
    	plt.ylabel('$\log_{10}\,X/X_\odot$')
    if ref == -2:
    	plt.ylabel('$\log_{10}\,X/X_{\mathrm{init}}$')
    if ref == -10:
    	plt.ylabel('$\log_{10}\,X$')
    if ref >= 0:
    	plt.ylabel('$\log_{10}\,X/X_{'+str(ref)+'}$')
    if labels == True:
    	for i in range(n_plot):
		if zrange != None and z_plot[i] >= zrange[0] and z_plot[i] <= zrange[-1]:
			dely = dely0
			extr = 0.0
    			if i%2 != 0:
				extr = extr0
			dely = dely + extr
			if i>0 and i<n_plot-2 and y_plot[i] < y_plot[i-1] and y_plot[i] < y_plot[i+1]:
        			dely = -0.2 - extr
			plt.text(z_plot[i],y_plot[i]+dely,el_name_plot[i],horizontalalignment='center',fontsize=label_fsize)
    plt.title(title)

    print (el_name_plot[78])

# The following function plots [X/Fe] decayed surface element ratios
# for a specified model from a specified mppnp H5_surf files
def plot_surf_xfe_abund(sol_ab,init_ab,dir,case,model,zrange=None,yrange=None,decay=True,shape='',msize=14,labels=True,label_fsize=12,title="",dely0=0.2,extr0=0.2,dilute=1.0):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    init_ab : string
       path to file with initial abundances that are used for scaling

    dir : string
       directory path to where H5_surf files are

    case : string
       H5_surf_case is the directory in the dir with surf files to plot

    model : integer
       model number to plot

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    decay : boolean
       if True then decayed surface abundance are plotted

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any              

    dely0 : float
       shift of element label in y

    extr0 : float
       additional shift of element label in y

    dilute : float
       dilution factor between 0. and 1.

    '''

# read in solar abundances in the path sol_ab
    utils.solar(sol_ab,1.)  
    sol_abu=utils.solar_elem_abund
    n_sol=len(sol_abu)

# read in initial abundances in the path init_ab
    utils.solar(init_ab,1.)  
    init_abu=utils.solar_elem_abund
    n_init=len(init_abu)

    if case == "":
    	s=nugridse.se(dir+"/H5_surf")
    else:
	s=nugridse.se(dir+"/H5_surf_"+case)

    if decay:

    	el_abu_hif=s.get(model,'elem_massf_decay')
    else:
    	el_abu_hif=s.get(model,'elem_massf')
    n_el_hif=len(el_abu_hif)
    z_hif=np.linspace(0,0,n_el_hif)
    y_hif=np.linspace(0,0,n_el_hif)
    el_name_hif=["  " for x in range(n_el_hif)]
    z_work=np.linspace(0,0,n_el_hif)
    y_work=np.linspace(0,0,n_el_hif)
    el_name_work=["  " for x in range(n_el_hif)]

    el_name_hif[0]='n'
    for i in range(n_el_hif):
    	z_hif[i]=float(i)       #  Z=i in mppnp surf data output
    	if (i>0):
        	el_name_hif[i]=utils.get_el_from_z(i)

    if dilute > 1.:
        dilute = 1.

    if dilute < 0.:
        dilute = 0.

# find delxfe = lg nFe/nFe_sol to change abundance ratios to [X/Fe]
    for i in range(n_el_hif):
        for k in range(n_sol):
            z_sol=k+1
            if float(z_sol)==z_hif[i] and z_sol == 1:
               delx = np.log10((dilute*el_abu_hif[i]+(1.-dilute)*init_abu[k])/sol_abu[k])
            if float(z_sol)==z_hif[i] and z_sol == 26:
               delxfe = np.log10((dilute*el_abu_hif[i]+(1.-dilute)*init_abu[k])/sol_abu[k])
    fehmod = delxfe - delx
    print "\n delxfe =",delxfe,", delx =",delx,", fehmod =",fehmod

# prepare plot arrays
    j=-1
    for i in range(n_el_hif):
    	for k in range(n_sol):
        	z_sol=k+1
       		if float(z_sol)==z_hif[i] and z_sol != 43 and z_sol != 61:
          		j=j+1
            		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if sol_abu[k]>1e-30 and el_abu_hif[i]>1e-30:
                		y_work[j]=np.log10((dilute*el_abu_hif[i]+(1.-dilute)*init_abu[k])/sol_abu[k])
    	n_plot=j+1

    z_plot=np.linspace(0,0,n_plot)
    y_plot=np.linspace(0,0,n_plot)
    el_name_plot=["  " for x in range(n_plot)]
    for i in range(n_plot):
    	z_plot[i]=z_work[i]
# delxfe changes abundance ratios to [X/Fe]
    	y_plot[i]=y_work[i] - delxfe
	el_name_plot[i]=el_name_work[i]

# plot [X/Fe] surface ratios mppnp hif abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':   
    	shape='k-'     
    plt.plot(z_plot,y_plot,shape,markersize=msize)
    #plt.plot(z_plot,y_plot,'k-')
    if zrange == None:
    	plt.xlim(0,85)
    else:
    	plt.xlim(zrange)
    if yrange == None:
    	plt.ylim(-5,5)
    else:
    	plt.ylim(yrange)
    plt.xlabel('charge number Z')
    plt.ylabel('[X/Fe]')
    if labels == True:
    	for i in range(n_plot):
		if zrange != None and z_plot[i] >= zrange[0] and z_plot[i] <= zrange[-1]:
			dely = dely0
			extr = 0.0
    			if i%2 != 0:
				extr = extr0
			dely = dely + extr
			if i>0 and i<n_plot-2 and y_plot[i] < y_plot[i-1] and y_plot[i] < y_plot[i+1]:
        			dely = -dely - extr
			plt.text(z_plot[i],y_plot[i]+dely,el_name_plot[i],horizontalalignment='center',fontsize=label_fsize)
    plt.title(title,fontsize=label_fsize)

    print (el_name_plot[78])

# The following function plots FIRST DERIVATIVES of the solar-scaled decayed surface element  
# abundances for a specified model from a specified mppnp H5_surf files
def plot_deriv_surf_abund(sol_ab,dir,case,model,ref=-1,zrange=None,yrange=None,shape='',msize=14,labels=True,title=""):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    ref : integer
       if ref > 0 then ref = model number to scale plotted abundances
       if ref = -1 then use sol_ab for the scaling
       if ref = -2 then still use sol_ab which are now in fact initial ones for the scaling 
       if ref = -10 then do not use any scaling

    dir : string
       directory path to where H5_surf files are

    case : string
       H5_surf_case is the directory in the dir with surf files to plot

    model : integer
       model number to plot

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any              
    '''

# read in solar abundances in the path sol_ab
    if ref == -1 or ref == -2:
    	utils.solar(sol_ab,1.)  
    	sol_abu=utils.solar_elem_abund
    	n_sol=len(sol_abu)

    if case == "":
    	s=nugridse.se(dir+"/H5_surf")
    else:
	s=nugridse.se(dir+"/H5_surf_"+case)

    el_abu_hif=s.get(model,'elem_massf_decay')
    n_el_hif=len(el_abu_hif)
    z_hif=np.linspace(0,0,n_el_hif)
    y_hif=np.linspace(0,0,n_el_hif)
    el_name_hif=["  " for x in range(n_el_hif)]
    z_work=np.linspace(0,0,n_el_hif)
    y_work=np.linspace(0,0,n_el_hif)
    el_name_work=["  " for x in range(n_el_hif)]

    el_name_hif[0]='n'
    for i in range(n_el_hif):
    	z_hif[i]=float(i)       #  Z=i in mppnp surf data output
    	if (i>0):
        	el_name_hif[i]=utils.get_el_from_z(i)

    if ref >= 0:
    	el_abu_ref=s.get(ref,'elem_massf_decay')
    	n_el_ref=len(el_abu_ref)
    	z_ref=np.linspace(0,0,n_el_ref)
    	y_ref=np.linspace(0,0,n_el_ref)
    	el_name_ref=["  " for x in range(n_el_ref)]

    	el_name_ref[0]='n'
    	for i in range(n_el_ref):
    		z_ref[i]=float(i)       #  Z=i in mppnp surf data output
    		if (i>0):
        		el_name_ref[i]=utils.get_el_from_z(i)

# prepare plot arrays
    if ref == -1 or ref == -2:
    	j=-1
    	for i in range(n_el_hif):
    		for k in range(n_sol):
        		z_sol=k+1
        		if float(z_sol)==z_hif[i] and z_sol != 43 and z_sol != 61:
            			j=j+1
            			z_work[j]=z_hif[i]
            			y_work[j]=-20.
				el_name_work[j]=el_name_hif[i]
            			if sol_abu[k]>1e-30 and el_abu_hif[i]>1e-30:
                			y_work[j]=np.log10(el_abu_hif[i]/sol_abu[k])
    	n_plot=j+1
    if ref >= 0:
    	j=-1
    	for i in range(n_el_hif):
		z_el=z_hif[i]
		if z_el != 0 and z_el != 43 and z_el != 61:
			j=j+1
        		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if el_abu_ref[i]>1e-30 and el_abu_hif[i]>1e-30:
               			y_work[j]=np.log10(el_abu_hif[i]/el_abu_ref[i])
	n_plot=j
    if ref == -10 :
    	j=-1
    	for i in range(n_el_hif):
		z_el=z_hif[i]
		if z_el != 0 and z_el != 43 and z_el != 61:
			j=j+1
        		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if el_abu_hif[i]>1e-30:
               			y_work[j]=np.log10(el_abu_hif[i])
	n_plot=j

    z_plot=np.linspace(0,0,n_plot)
    y_plot=np.linspace(0,0,n_plot)
    el_name_plot=["  " for x in range(n_plot)]
#   for i in range(n_plot):
    for i in range(n_plot-1):
    	z_plot[i]=z_work[i]
    	y_plot[i]=y_work[i+1] - y_work[i]
	el_name_plot[i]=el_name_work[i]

# plot first derivatives of the solar scaled surface mppnp hif abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':   
    	shape='bo'     
    plt.plot(z_plot,y_plot,shape,markersize=msize)
    plt.plot(z_plot,y_plot,'k-')
    if zrange == None:
    	plt.xlim(0,85)
    else:
    	plt.xlim(zrange)
    if yrange == None:
    	plt.ylim(-5,5)
    else:
    	plt.ylim(yrange)
    plt.xlabel('charge number Z')
    if ref == -1:
    	plt.ylabel('$\Delta \log_{10}\,X/X_\odot$')
    if ref == -2:
    	plt.ylabel('$\Delta \log_{10}\,X/X_{\mathrm{init}}$')
    if ref == -10:
    	plt.ylabel('$\Delta \log_{10}\,X$')
    if ref >= 0:
    	plt.ylabel('$\Delta \log_{10}\,X/X_{'+str(ref)+'}$')
    if labels == True:
    	for i in range(n_plot):
    		if i%2 != 0:
        		plt.text(z_plot[i],y_plot[i]+0.2,el_name_plot[i])
    		else:
#       		plt.text(z_plot[i],y_plot[i]-0.2,el_name_plot[i])
        		plt.text(z_plot[i],y_plot[i]+0.2,el_name_plot[i])
    plt.title(title)

    print (el_name_plot[78])

# The following function plots solar-scaled decayed element  
# abundances for a specified cycle (time step) from a specified directory with results of PPN simulations
def plot_ppn_abund(sol_ab,ppn_dir,cycle,ref=-1,zrange=None,yrange=None,clr='',mrk='o',mrkclr='k',shape='',msize=14,labels=True,label_fsize=12,title=""):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    ref : integer
       if ref > 0 then ref = cycle number to scale plotted abundances
       if ref = -1 then use sol_ab for the scaling
       if ref = -2 then still use sol_ab which are now in fact initial ones for the scaling 
       if ref = -10 then do not use any scaling

    ppn_dir : string
       directory path to where results of PPN simulations are

    cycle : integer
       cycle number to plot

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any      
    '''

# read in solar abundances in the path sol_ab
    if ref == -1 or ref == -2:
    	utils.solar(sol_ab,1.)
    	sol_abu=utils.solar_elem_abund
    	n_sol=len(sol_abu)

    pb=ppn.abu_vector(ppn_dir)  # location of ppn simulation results

# this is a fragment of elemental_abund function from data_plot.py that prepares arrays for the plotting
    ztest=[0,200]
    title_items=None
    pb.get(cycle,decayed=True)
    z_el=np.unique(pb.z_iso_to_plot)
    zmin_ind=min(np.where(z_el>=ztest[0])[0])
    zmax_ind=max(np.where(z_el<=ztest[1])[0])

# extract some elemental quantities:
    a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
    for z in z_el[zmin_ind:zmax_ind]:
    	el=pb.el_iso_to_plot[np.where(pb.z_iso_to_plot==z)[0].tolist()[0]]
    	X_el=pb.abunds[np.where(pb.el_iso_to_plot==el)[0].tolist()].sum()
    	a_el.append(pb.a_iso_to_plot[np.where(pb.z_iso_to_plot==z)[0].tolist()[0]])
    	el_abu.append(X_el)
    	el_name.append(el)
    	el_abu_hash[el]=X_el

    n_el=len(z_el)-1  # number of stable (decayed) elements to plot
    z_el[n_el]  # the maximum Z in the results of ppn simulations

# extract reference elemental abundances:
    if ref >= 0:
    	pb_ref=ppn.abu_vector(ppn_dir)  # location of ppn simulation results
    	pb_ref.get(ref,decayed=True)
    	el_abu_ref=[]
    	for z in z_el[zmin_ind:zmax_ind]:
    		el=pb_ref.el_iso_to_plot[np.where(pb_ref.z_iso_to_plot==z)[0].tolist()[0]]
    		X_el=pb_ref.abunds[np.where(pb_ref.el_iso_to_plot==el)[0].tolist()].sum()
    		el_abu_ref.append(X_el)

    z=np.linspace(0,0,n_el)
    y=np.linspace(0,0,n_el)

    if ref == -1 or ref == -2:
    	for i in range(n_el):
    		for k in range(n_sol):
        		z_sol=k+1
        		if z_sol == z_el[i]:
            			z[i]=z_el[i]
            			y[i]=-20.
            			if sol_abu[k]>1e-30:
                			y[i]=np.log10(el_abu[i]/sol_abu[k])
    if ref >= 0:
        for i in range(n_el):
            	z[i]=z_el[i]
        	y[i]=np.log10(el_abu[i]/el_abu_ref[i])
    if ref == -10:
        for i in range(n_el):
            	z[i]=z_el[i]
        	y[i]=np.log10(el_abu[i])

# plot scaled PPN abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':
        plt.plot(z,y,color=clr,marker=mrk,markerfacecolor=mrkclr,markersize=msize)
    else:
        plt.plot(z,y,shape,markersize=msize)
    #plt.plot(z_plot,y_plot,'k-')
    if zrange == None:
        plt.xlim(0,85)
    else:
        plt.xlim(zrange)
    if yrange == None:
        plt.ylim(-5,5)
    else:
        plt.ylim(yrange)
    plt.xlabel('charge number Z')
    if ref == -1:
    	plt.ylabel('$\log_{10}\,X/X_\odot$')
    if ref == -2:
    	plt.ylabel('$\log_{10}\,X/X_{\mathrm{init}}$')
    if ref == -10:
    	plt.ylabel('$\log_{10}\,X$')
    if ref >= 0:
    	plt.ylabel('$\log_{10}\,X/X_{'+str(ref)+'}$')
    if labels == True:
        for i in range(n_el):
		if zrange != None and z[i] >= zrange[0] and z[i] <= zrange[-1]:
			dely = 0.2
			extr = 0.0
                	if i%2 != 0:
				extr = 0.2
			dely = dely + extr
			if i>1 and i<n_el-2 and y[i] < y[i-1] and y[i] < y[i+1]:
        			dely = -0.2 - extr
                        #plt.text(z[i],y[i]+dely,el_name[i],fontsize=label_fsize)  # i+1 in the el name to skip neutron
                        plt.text(z[i],y[i]+dely,el_name[i],ha='center',fontsize=label_fsize)  # i+1 in the el name to skip neutron
    plt.title(title)

# The following function plots [X/Fe] decayed element ratios  
# for a specified cycle (time step) from a specified directory with results of PPN simulations
def plot_ppn_xfe_abund(sol_ab,ppn_dir,cycle,zrange=None,yrange=None,shape='',msize=14,labels=True,label_fsize=12,title="",dely0=0.2,extr0=0.2,dilute=1.0,feh=0.0):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    ppn_dir : string
       directory path to where results of PPN simulations are

    cycle : integer
       cycle number to plot

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any      

    dely0 : float
       shift of element label in y

    extr0 : float
       additional shift of element label in y

    dilute : float
       decimal logarithm of the dilution factor

    feh : float
       [Fe/H] ratio                      
    '''

# read in solar abundances in the path sol_ab
    utils.solar(sol_ab,1.)
    sol_abu=utils.solar_elem_abund
    n_sol=len(sol_abu)

    pb=ppn.abu_vector(ppn_dir)  # location of ppn simulation results

# this is a fragment of elemental_abund function from data_plot.py that prepares arrays for the plotting
    ztest=[0,200]
    title_items=None
    pb.get(cycle,decayed=True)
    z_el=np.unique(pb.z_iso_to_plot)
    zmin_ind=min(np.where(z_el>=ztest[0])[0])
    zmax_ind=max(np.where(z_el<=ztest[1])[0])

# extract some elemental quantities:
    a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
    for z in z_el[zmin_ind:zmax_ind]:
    	el=pb.el_iso_to_plot[np.where(pb.z_iso_to_plot==z)[0].tolist()[0]]
    	X_el=pb.abunds[np.where(pb.el_iso_to_plot==el)[0].tolist()].sum()
    	a_el.append(pb.a_iso_to_plot[np.where(pb.z_iso_to_plot==z)[0].tolist()[0]])
    	el_abu.append(X_el)
    	el_name.append(el)
    	el_abu_hash[el]=X_el

    n_el=len(z_el)-1  # number of stable (decayed) elements to plot
    z_el[n_el]  # the maximum Z in the results of ppn simulations

    z=np.linspace(0,0,n_el)
    y=np.linspace(0,0,n_el)

# it is assumed that X = X_sol
# take into account a reduced initial [Fe/H] and dilution
    delxfe = feh - dilute

    for i in range(n_el):
    	for k in range(n_sol):
       		z_sol=k+1
       		if z_sol == z_el[i]:
         		z[i]=z_el[i]
            		y[i]=-20.
            		if sol_abu[k]>1e-30:
                		y[i]=np.log10(el_abu[i]/sol_abu[k]) - delxfe

# plot scaled PPN abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':
        shape='k-'
    plt.plot(z,y,shape,markersize=msize)
    #plt.plot(z,y,'k-')
    if zrange == None:
        plt.xlim(0,85)
    else:
        plt.xlim(zrange)
    if yrange == None:
        plt.ylim(-5,5)
    else:
        plt.ylim(yrange)
    plt.xlabel('charge number Z')
    plt.ylabel('[X/Fe]')
    if labels == True:
        for i in range(n_el):
		if zrange != None and z[i] >= zrange[0] and z[i] <= zrange[-1]:
			dely = dely0
			extr = 0.0
                	if i%2 != 0:
				extr = extr0
			dely = dely + extr
			if i>1 and i<n_el-2 and y[i] < y[i-1] and y[i] < y[i+1]:
        			dely = -dely - extr
                        plt.text(z[i],y[i]+dely,el_name[i],fontsize=label_fsize)  # i+1 in the el name to skip neutron
    plt.title(title)

def plot_isoabund(p,cycle, ifig=151, stable = False, decayed = False, \
                  amass_range = None, ylim = [0,0], ref = -1, title_format = 1, 
                  img_title = "", save_fig=False):
    """
    plots the abundance of all the chemical species for a range of cycles

    input:
    p: abu_vector instance
    cycle: cycle for which to plot
    stable: a boolean of whether to filter out the unstables.
        Defaults to False
    decayed: [True/]False to plot decayed distributions (True), or life 
        distribution
    amass_range: a 1x2 array containing the lower and upper Atomic
        mass range. optional. if None plot entire available
        atomic mass range
    ylim: A 1x2 array containing the lower and upper Y limits.
        Defaults to [0,0], in which case ylim will be determined automatically
    ref: reference cycle. If it is not -1 (default), this method will
        plot the abundences of cycle devided by the cycle of the same
        instance given in the ref variable. If ref is a list it will 
        be interpreted to have two elements: ref=['dir/of/ref/run',cycle]
        which uses a refernece cycle from another run. If any
        abundence in the reference cycle is zero, it will replace it with
        1e-99. The default is -1, it will do nothing.
    title_format: changes the formating of the title, defaults to 1
    img_title: If it is not "" img_title sets title of the immage, if
        it is "" the title is generated automatically based on the
        title_format, defaults to ""
    save_fig: save image to file
    """

    mp=p.get('mod')
    form_str='%6.1F'
    form_str1='%4.3F'
    #ylimits=(-11,-3)#var
    i=0
    T9=p.get('t9',fname=cycle)
    Rho=p.get('rho',fname=cycle)
    mod=p.get('mod',fname=cycle)
    time= p.get('agej',fname=cycle)*utils.constants.one_year
    plt.close(ifig); plt.figure(ifig)

    p.iso_abund(cycle,amass_range=amass_range,decayed=decayed,stable=stable,ylim=ylim,ref=ref,show=False)
    if img_title == "":
        if title_format == 1:
            plt.title(str(mod)+' t='+form_str%time+'s $T_9$='+form_str1%T9+' $\\rho$='+str(Rho))
        elif title_format == 2:
            plt.title('model='+str(mod)+' T9='+str(T9)+' Rho='+str(Rho))
    else:
        plt.title(img_title)
    if save_fig:
        plt.savefig('iso_abund_'+data_plot._padding_model_number(cycle,mp[len(mp)-1])+'.png')

def plot_abu_evolution(x,ifig=147,specs = ['PROT','C  12','C  13','N  13','N  14','O  16','KR 89','I 135'],
                       y_lim = (-6,0.2), legend_loc = 0, x_axis_rev = False, log_time = False,
                       y_axis_offset = 0, time_in_min = True, extended_label = "",save_fig=False):
    """
    Plots the abundance of species over time

    input:
    x    : xtime instance
    specs: The list of species to track from the x-time file, defaults
        to ['PROT','HE  4','C  12','N  14','O  16','NE 20','MG 24','SI 28']
    y_lim: The range of the y axis, defaults to (-5, 0.2)
    legend_loc: The location of the legend on the plot, defaults to 4
    x_axis_rev: turns off/on the revircing of the x axis, defaults
        to False
    log_time: turns off/on logrithmic time, defaults to False
    y_axis_offset: applyes a constant offset to the y axis,
        defaults to 0
    time_in_min: turns off/on the display of time in min vs years,
        defaults to True
    extended_label: A string with extra information to append on to
        the end of the label, defaults to ""    
    save_fig      : default to False, save figure to 'abu_evolution.png' 
    """

    abus=[]
    for spec in specs:
        abu=x.get(spec)
        abus.append(abu)

    if time_in_min:
        time=x.get('time')*(3.1558149984e7/60.)
    else:
       time=x.get('time')
    plt.close(ifig);plt.figure(ifig)
    for i in range(len(specs)):
        shape,mevery=utils.linestyle(i,50,10)
        if log_time:
            plt.plot(np.log10(time),np.log10(abus[i] + y_axis_offset),shape,markevery=mevery,label=specs[i])
        else:
            plt.plot(time,np.log10(abus[i] + y_axis_offset),shape,markevery=mevery,label=specs[i])

    if x_axis_rev:
        utils.xlimrev()
    plt.ylim(*y_lim)
    plt.legend(loc=legend_loc)
    x_label = ''

    if log_time:
        x_label += '$\log t '
    else:
        x_label += '$t / '
        if time_in_min:
            x_label += '\mathrm{min}'
        else:
            x_label += '\mathrm{yrs}'

    x_label += extended_label
    x_label += "$"
    plt.xlabel(x_label)
    plt.ylabel('$\log X \mathrm{[mass fraction]}$')
    if save_fig:
        plt.savefig('abu_evolution.png')


def t_tau(x):
    '''plot neutron exposure as a function of time
    
    Parameter : instance
     x  a ppn.xtime instance

    '''

    xn=x.get("NEUT")
    age=x.get("time")
    rho=x.get("rho")
    v_therm=3.e8
    oneyear=utils.constants.one_year
    #Nn=1.e-19*ast.avogadro_constant*1000.
    #Nn/1.e7

    Nn=xn*ast.avogadro_constant*rho
    age_s=age*oneyear
    tau=integrate.cumtrapz(Nn*v_therm,age_s)
    tau=np.log10(tau*1.e-27)

    ifig=5;plt.close(ifig);plt.figure(ifig)
    plt.plot(np.log10(age[1:]),tau)
    plt.xlabel('$\log_{10}\, t\ (\mathrm{yr})$')
    plt.ylabel('$\tau\ (\mathrm{[mbarn^{-1}]})$')


    ifig=4;plt.close(ifig);plt.figure(ifig)
    plt.plot(x.get('cycle')[1:],tau)
    plt.xlabel('model number')
    plt.ylabel('$\\tau$ $\mathrm{[mbarn^{-1}]}$')

    ifig=3;plt.close(ifig);plt.figure(ifig)
    plt.plot(x.get('cycle'),np.log10(Nn))
    plt.xlabel('model number')
    plt.ylabel('neutron density $\\log N_\mathrm{n} \mathrm{cm}^{-3}$')
    plt.ylim(4,17)

    ifig=2;plt.close(ifig);plt.figure(ifig)
    plt.plot(np.log10(age),np.log10(Nn))
    plt.ylabel('neutron density $\\log N_\mathrm{n} \mathrm{cm}^{-3}$')
    plt.xlabel('$\\log t/\mathrm{yr}$')
    plt.ylim(4,17)

def elemental(p,cycles_input,zrange=(10,30),yrange=(-1,3),\
              save_fig=False,show=True):
    '''Make series of elemental overabundance plots

    p              abu_vector instance
    cycles triple: start_cycle, stop_cycle, stride
    zrange         range of atomic numbers to plot
    yrange         yrange (log scale, [X/Fe])
    show           show image on screen (True, default) or save images in png file

    Note:
    This methods needs to become part of data_plot and be called directly via 
    ppn.py/nugrdise.py!
    '''

    print ("This is elemental in local_utils.")

    stride = cycles_input[2]
    cycles=range(cycles_input[0],cycles_input[1],stride)
    form_str='%4.1F'
    if min(cycles) is not 0:  # making sure the first one is there 
        cycles = [0] + cycles
    max_num = max(cycles) # holding the initial abundance
    case_name = p.sldir.split('/')[-2] # name of case, hack

    for i in cycles:
        p.get(i,decayed=True) # this adds decayed 'abunds' elemental array to instance
        # list of unique charges/elements: 
        z_el=np.unique(p.z_iso_to_plot)
        zmin=min(np.where(z_el>zrange[0])[0])
        zmax=max(np.where(z_el<zrange[1])[0])-1

        # extract some elemental quantities:
        a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
        for z in z_el[zmin:zmax]:
            el=p.el_iso_to_plot[np.where(p.z_iso_to_plot==z)[0].tolist()[0]]
            X_el=p.abunds[np.where(p.el_iso_to_plot==el)[0].tolist()].sum()
            a_el.append(p.a_iso_to_plot[np.where(p.z_iso_to_plot==z)[0].tolist()[0]])
            el_abu.append(X_el)
            el_name.append(el)
            el_abu_hash[el]=X_el

        if i == 0:
            el_abu0 = np.copy(el_abu)
            continue
        el_abu /= el_abu0 # normalize to initial abundance
                          # ! here we need to add the option to divide by
                          # ! any solar abundance distribution from USEEPP

        # prepare canvas
        istr = str(i).zfill(len(str(max_num)))
        if show:
            plt.figure(i)
        else:
            plt.close(1);plt.figure(1)

        # plot an elemental abundance distribution with labels:
        plt.plot(z_el[zmin:zmax],np.log10(el_abu),'o-')
        j=0        # add labels
        for z in z_el[zmin:zmax]:
            plt.text(z-0.5,np.log10(el_abu[j])+0.05,el_name[j])
            j += 1

        # neutron density:
        species='N-1'
        xn=p.get('ABUNDANCE_MF',i)[int(np.where(p.get('ISOTP',i)==species)[0].var())]
        rho = p.get('rho',i)
        v_therm=3.e8
        Nn=xn*ast.avogadro_constant*rho
        oneyear=utils.constants.one_year
        age_min=p.get('agej',i)*oneyear/60.
        plt.title('run: '+case_name+', cycle: '+istr+\
                  ', log Nn='+form_str%np.log10(Nn)+', t/min='\
                  +form_str%age_min )
        plt.ylim(yrange)
        plt.xlim(zrange) 
        plt.xlabel('charge number Z')
        plt.ylabel('$log El / El_\odot$')   
        plt.legend(loc=0)
        if save_fig: 
            plt.savefig('overabund_'+istr+'.png')

class one_zone():
    '''
    read case summary of single zone calcualtions and provide initialization 
    routine for cases
    '''
    def __init__(self,data_dir, case_summary):
        '''Initialize

        data_dir        where the one-zone calcualtions and the summary file is
        case_summary    name of summary file, format according to "1-zone models" notebook
        '''

        self.data_dir=data_dir

        ff=open(data_dir+case_summary)
        all_lines=ff.readlines()
        self.cases=[];self.dir_names=[];self.initial_H=[];self.initial_Z=[]
        self.T6=[];self.rho=[];self.comments=[]
        for line in all_lines:
            case,this_dir,Hini,Zini,this_T6,this_rho,this_comment = line.split('|')
            self.cases.append(case.strip()); self.dir_names.append(this_dir.strip())
            self.initial_H.append(float(Hini))
            self.initial_Z.append(float(Zini)); self.T6.append(float(this_T6))
            self.rho.append(float(this_rho)); self.comments.append(this_comment.strip())
        ff.close() 

    def initialize_ppn_case(self,case,data_type='xtime'):
        '''create ppn instance of data_type for case

        case string: one of the available cases in array cases
        data_type    either xtime or abu_vector
        '''
        import ppn
        if data_type is 'xtime':
            return ppn.xtime(self.data_dir+self.dir_names[self.cases.index(case)])
        elif data_type is 'abu_vector':
            return ppn.abu_vector(self.data_dir+self.dir_names[self.cases.index(case)])

# The following function plots solar-scaled decayed surface element  
# abundances for a specified model from a specified mppnp H5_surf files
# this is a modification of plot_surf_abund that assumes that 90Sr has not decayed and, as a result,
# the final Sr is higher by +0.42dex and final Zr is lower by -0.16dex
def plot_surf_abund_corr(sol_ab,dir,case,model,ref=-1,zrange=None,yrange=None,corrSr=0.42,corrZr=-0.16,shape='',msize=14,labels=True,title=""):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    ref : integer
       if ref > 0 then ref = model number to scale plotted abundances
       if ref = -1 then use sol_ab for the scaling
       if ref = -2 then still use sol_ab which are now in fact initial ones for the scaling 
       if ref = -10 then do not use any scaling

    dir : string
       directory path to where H5_surf files are

    case : string
       H5_surf_case is the directory in the dir with surf files to plot

    model : integer
       model number to plot

    zrange : tuple, float, optional 
       range in Z for plot

    yrange : tuple, float, optional 
       range in Y for plot

    corrSr, corrZr : float   
	correction (in dex) to the Sr and Zr abundances

    shape : string
       plotting symbol and color

    labels : boolean
       if True then element labels are plotted

    title : string
       title to plot if any              
    '''

# read in solar abundances in the path sol_ab
    if ref == -1 or ref == -2:
    	utils.solar(sol_ab,1.)  
    	sol_abu=utils.solar_elem_abund
    	n_sol=len(sol_abu)

    if case == "":
    	s=nugridse.se(dir+"/H5_surf")
    else:
	s=nugridse.se(dir+"/H5_surf_"+case)

    el_abu_hif=s.get(model,'elem_massf_decay')
    n_el_hif=len(el_abu_hif)
    z_hif=np.linspace(0,0,n_el_hif)
    y_hif=np.linspace(0,0,n_el_hif)
    el_name_hif=["  " for x in range(n_el_hif)]
    z_work=np.linspace(0,0,n_el_hif)
    y_work=np.linspace(0,0,n_el_hif)
    el_name_work=["  " for x in range(n_el_hif)]

    el_name_hif[0]='n'
    for i in range(n_el_hif):
    	z_hif[i]=float(i)       #  Z=i in mppnp surf data output
    	if (i>0):
        	el_name_hif[i]=utils.get_el_from_z(i)

    if ref >= 0:
    	el_abu_ref=s.get(ref,'elem_massf_decay')
    	n_el_ref=len(el_abu_ref)
    	z_ref=np.linspace(0,0,n_el_ref)
    	y_ref=np.linspace(0,0,n_el_ref)
    	el_name_ref=["  " for x in range(n_el_ref)]

    	el_name_ref[0]='n'
    	for i in range(n_el_ref):
    		z_ref[i]=float(i)       #  Z=i in mppnp surf data output
    		if (i>0):
        		el_name_ref[i]=utils.get_el_from_z(i)

# prepare plot arrays
    if ref == -1 or ref == -2:
    	j=-1
    	for i in range(n_el_hif):
    		for k in range(n_sol):
        		z_sol=k+1
        		if float(z_sol)==z_hif[i] and z_sol != 43 and z_sol != 61:
            			j=j+1
            			z_work[j]=z_hif[i]
            			y_work[j]=-20.
				el_name_work[j]=el_name_hif[i]
            			if sol_abu[k]>1e-30 and el_abu_hif[i]>1e-30:
                			y_work[j]=np.log10(el_abu_hif[i]/sol_abu[k])
    	n_plot=j+1
    if ref >= 0:
    	j=-1
    	for i in range(n_el_hif):
		z_el=z_hif[i]
		if z_el != 0 and z_el != 43 and z_el != 61:
			j=j+1
        		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if el_abu_ref[i]>1e-30 and el_abu_hif[i]>1e-30:
               			y_work[j]=np.log10(el_abu_hif[i]/el_abu_ref[i])
	n_plot=j
    if ref == -10 :
    	j=-1
    	for i in range(n_el_hif):
		z_el=z_hif[i]
		if z_el != 0 and z_el != 43 and z_el != 61:
			j=j+1
        		z_work[j]=z_hif[i]
            		y_work[j]=-20.
			el_name_work[j]=el_name_hif[i]
            		if el_abu_hif[i]>1e-30:
               			y_work[j]=np.log10(el_abu_hif[i])
	n_plot=j

    z_plot=np.linspace(0,0,n_plot)
    y_plot=np.linspace(0,0,n_plot)
    el_name_plot=["  " for x in range(n_plot)]
    for i in range(n_plot):
    	z_plot[i]=z_work[i]
    	y_plot[i]=y_work[i]
	if z_plot[i] == 38:
		y_plot[i] = y_plot[i] + corrSr
	if z_plot[i] == 40:
		y_plot[i] = y_plot[i] + corrZr
	el_name_plot[i]=el_name_work[i]

# plot solar scaled surface mppnp hif abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':   
    	shape='bo'     
    plt.plot(z_plot,y_plot,shape,markersize=msize)
    plt.plot(z_plot,y_plot,'k-')
    if zrange == None:
    	plt.xlim(0,85)
    else:
    	plt.xlim(zrange)
    if yrange == None:
    	plt.ylim(-5,5)
    else:
    	plt.ylim(yrange)
    plt.xlabel('charge number Z')
    if ref == -1:
    	plt.ylabel('$\log_{10}\,X/X_\odot$')
    if ref == -2:
    	plt.ylabel('$\log_{10}\,X/X_{\mathrm{init}}$')
    if ref == -10:
    	plt.ylabel('$\log_{10}\,X$')
    if ref >= 0:
    	plt.ylabel('$\log_{10}\,X/X_{'+str(ref)+'}$')
    if labels == True:
    	for i in range(n_plot):
    		if i%2 != 0:
        		plt.text(z_plot[i],y_plot[i]+0.2,el_name_plot[i])
    		else:
#       		plt.text(z_plot[i],y_plot[i]-0.2,el_name_plot[i])
        		plt.text(z_plot[i],y_plot[i]+0.2,el_name_plot[i])
    plt.title(title)

    print (el_name_plot[78])

# The following function plots one abundance ratio vs. another, e.g. [Rb/Sr] vs. [Y/Sr]
# for a specified cycle (time step) from a specified directory with results of PPN simulations
def plot_ppn_rel_abund(sol_ab,ppn_dir,cycle,z1plot=37,z2plot=38,z3plot=39,xxrange=None,yyrange=None,shape='',msize=14,title=""):
    '''
    Parameters 
    ----------
    sol_ab : string
       path to file with solar abundances that are used for scaling

    ppn_dir : string
       directory path to where results of PPN simulations are

    cycle : integer
       cycle number to plot

    xxrange : tuple, float, optional 
       range in X for plot

    yyrange : tuple, float, optional 
       range in Y for plot

    shape : string
       plotting symbol and color

    title : string
       title to plot if any      
    '''

# read in solar abundances in the path sol_ab
    utils.solar(sol_ab,1.)
    sol_abu=utils.solar_elem_abund
    n_sol=len(sol_abu)

    pb=ppn.abu_vector(ppn_dir)  # location of ppn simulation results

# this is a fragment of elemental_abund function from data_plot.py that prepares arrays for the plotting
    ztest=[0,200]
    title_items=None
    pb.get(cycle,decayed=True)
    z_el=np.unique(pb.z_iso_to_plot)
    zmin_ind=min(np.where(z_el>=ztest[0])[0])
    zmax_ind=max(np.where(z_el<=ztest[1])[0])

# extract some elemental quantities:
    a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
    for z in z_el[zmin_ind:zmax_ind]:
    	el=pb.el_iso_to_plot[np.where(pb.z_iso_to_plot==z)[0].tolist()[0]]
    	X_el=pb.abunds[np.where(pb.el_iso_to_plot==el)[0].tolist()].sum()
    	a_el.append(pb.a_iso_to_plot[np.where(pb.z_iso_to_plot==z)[0].tolist()[0]])
    	el_abu.append(X_el)
    	el_name.append(el)
    	el_abu_hash[el]=X_el

    n_el=len(z_el)-1  # number of stable (decayed) elements to plot
    z_el[n_el]  # the maximum Z in the results of ppn simulations

    z=np.linspace(0,0,n_el)
    y=np.linspace(0,0,n_el)

    for i in range(n_el):
    	for k in range(n_sol):
       		z_sol=k+1
       		if z_sol == z_el[i]:
         		z[i]=z_el[i]
            		y[i]=-20.
            		if sol_abu[k]>1e-30:
                		y[i]=np.log10(el_abu[i]/sol_abu[k])

    xx=np.linspace(0,0,1)
    yy=np.linspace(0,0,1)

    for i in range(n_el):
	if z[i] == z1plot:
                i1plot = i
		y1plot = y[i]
	if z[i] == z2plot:
                i2plot = i
		y2plot = y[i]
	if z[i] == z3plot:
                i3plot = i
		y3plot = y[i]
    
    xx[0] = y3plot - y2plot
    yy[0] = y1plot - y2plot

# plot scaled PPN abundances
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 16}
    plt.rc('font', **font)

    if shape == '':
        shape='bo'
#   plt.plot(xx,yy,shape,markersize=msize)
    plt.plot(xx,yy,marker=r"$ {} $".format(shape),markersize=msize)
    if xxrange == None:
        plt.xlim(-1.5,1.5)
    else:
        plt.xlim(xxrange)
    if yyrange == None:
        plt.ylim(-1,2)
    else:
        plt.ylim(yyrange)
    plt.xlabel('['+el_name[i3plot]+'/'+el_name[i2plot]+']')
    plt.ylabel('['+el_name[i1plot]+'/'+el_name[i2plot]+']')
    plt.title(title)

# the following function plots one relative abundance vs. another, e.g. [Rb/Sr] vs [Y/Sr]
def plot_asplund99_sakurai_rel_abund(data_dir,z1plot=37,z2plot=38,z3plot=39,xxrange=None,yyrange=None,msize=12):
    '''
    read in Asplund's 1999 Sakurai's element abundances (12.0 for H)

    Parameters 
    ----------
    data_dir : string
       directory path to where the file Asplund1999_Sakurai.txt is

    zrange : tuple, float, optional 
       range in Z for plot
    '''
    
    lines=open(data_dir+'Asplund1999_Sakurai.txt','r').readlines()
    nelsak=np.size(lines)
    zsak=np.linspace(0,0,nelsak)
    elsol=np.linspace(0,0,nelsak)
    elname=["  " for x in range(nelsak)]
    elsak=np.linspace(0,0,nelsak)
    i=0
    for line in lines:
	parts=line.split()
 	zsak[i]=np.array(float(parts[0]))
 	elsol[i]=np.array(float(parts[1]))
	elname[i]=parts[2]
 	elsak[i]=np.array(float(parts[3]))
	i += 1
    for i in range(nelsak):
        if zsak[i] == 26:
            isakfe = i
    y=np.linspace(0,0,nelsak)
    feh=elsak[isakfe]-elsol[isakfe]
    for i in range(nelsak):
#	y[i]=elsak[i]-elsol[i]-feh
	y[i]=elsak[i]-elsol[i]
    xerr=np.linspace(0,0,nelsak)
    yerr=np.linspace(0,0,nelsak)
    for i in range(nelsak):
   	yerr[i]=0.3 

    xx=np.linspace(0,0,1)
    yy=np.linspace(0,0,1)
    exx=np.linspace(0,0,1)
    eyy=np.linspace(0,0,1)

    for i in range(nelsak):
	if zsak[i] == z1plot:
                i1plot = i
		y1plot = y[i]
	if zsak[i] == z2plot:
                i2plot = i
		y2plot = y[i]
	if zsak[i] == z3plot:
                i3plot = i
		y3plot = y[i]
    
    xx[0] = y3plot - y2plot
    yy[0] = y1plot - y2plot
    exx[0] = 0.6
    eyy[0] = 0.6

    plt.errorbar(xx,yy,exx,eyy,'ko',markersize=msize)

    if xxrange == None:
        plt.xlim(-1.5,1.5)
    else:
        plt.xlim(xxrange)
    if yyrange == None:
        plt.ylim(-1,2)
    else:
        plt.ylim(yyrange)
    plt.xlabel('['+elname[i3plot]+'/'+elname[i2plot]+']')
    plt.ylabel('['+elname[i1plot]+'/'+elname[i2plot]+']')
