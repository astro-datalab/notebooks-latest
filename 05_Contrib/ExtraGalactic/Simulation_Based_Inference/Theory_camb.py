import camb
import healpy as hp
import numpy as np
from astropy.io import fits
from tabulate import tabulate
from camb import model, initialpower
import frogress

class cosmo(object):
    def __init__(self,H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.0, As=2e-9, ns=0.965, r=0, num_massive_neutrinos=0):
        self.H0    = H0
        self.h     = H0/100.
        self.ombh2 = ombh2
        self.omch2 = omch2
        self.omc   = (omch2)/self.h**2
        self.omb   = (ombh2)/self.h**2
        self.mnu   = mnu
        self.omk   = omk
        self.tau   = tau
        self.As    = As
        self.ns    = ns
        self.r     = r
        self.num_massive_neutrinos =num_massive_neutrinos

class theory(object):
    """ This class is used to calculate the theory curves with limber"""
    def __init__(self, sigma_8 = None, cosmo=None,lmax=4000,nz=140,fast=False,feedback=1,halofit_version='mead',cambini=None, chistar = None, max_zs = 12.):

        self.lmax  = lmax
        self.nz    = nz
        self.h     = cosmo.h
        self.ombh2 = cosmo.ombh2
        self.omch2 = cosmo.omch2

        camb.set_feedback_level(level=feedback)

        if cosmo!=None:
            if fast==True:
                print("Using fast settings. Results my not be accurate")
                AccuracyBoost=1; max_eta_k=7000; lens_potential_accuracy=1; kmax=5
            else:
                AccuracyBoost=2; max_eta_k=50000; lens_potential_accuracy=4; kmax=500
            self.pars = camb.CAMBparams()
            self.pars.set_cosmology(H0=cosmo.H0, ombh2=cosmo.ombh2, omch2=cosmo.omch2, mnu=cosmo.mnu, omk=cosmo.omk, tau=cosmo.tau, num_massive_neutrinos=cosmo.num_massive_neutrinos)
            self.pars.InitPower.set_params(As=cosmo.As, ns=cosmo.ns, r=cosmo.r)
            self.pars.set_for_lmax(lmax=self.lmax,lens_potential_accuracy=lens_potential_accuracy, max_eta_k=max_eta_k)
            self.pars.set_accuracy(AccuracyBoost=AccuracyBoost,lSampleBoost=1, lAccuracyBoost=1)
            self.pars.NonLinearModel.set_params(halofit_version=halofit_version)
            
            if sigma_8 != None:
                self.pars.set_matter_power(redshifts = [0.], kmax = 2.0)
                results_ = camb.get_results(self.pars)
                r0 = (sigma_8**2/ results_.get_sigma8()**2)
                self.pars.InitPower.set_params(As=cosmo.As*r0, ns=cosmo.ns, r=cosmo.r)
                self.pars.set_for_lmax(lmax=self.lmax,lens_potential_accuracy=lens_potential_accuracy, max_eta_k=max_eta_k)
                self.pars.set_accuracy(AccuracyBoost=AccuracyBoost,lSampleBoost=1, lAccuracyBoost=1)
                self.pars.NonLinearModel.set_params(halofit_version=halofit_version)
                self.pars.set_matter_power(redshifts = [0.], kmax = 2.0)
                results_ = camb.get_results(self.pars)
                
                    
            self.results = camb.get_results(self.pars)
            
            
            self.h     = cosmo.h
            self.ombh2 = cosmo.ombh2
            self.omch2 = cosmo.omch2
            camb.set_feedback_level(level=feedback)

        if cambini!=None:
            self.pars    = camb.read_ini(cambini)
            self.h       = self.pars.H0/100.
            self.ombh2   = self.pars.ombh2
            self.omch2   = self.pars.omch2
            self.results = camb.get_results(self.pars)

        if ((cambini!=None) & (cosmo!=None)):
            sys.exit("need to choose cosmology or feed camb .ini file")

        self.omm   = (self.pars.omch2+self.pars.ombh2)/(self.h)**2


        #-----------------------------------------------------------------------
        # constants
        self.c       = 299792458.
        self.fc      = 1.5*(1000*100)**2*(self.ombh2/self.h**2+self.omch2/self.h**2)/self.c/self.c
        
        #-----------------------------------------------------------------------
        # lin version
        self.chistar = (self.results.conformal_time(0)- self.results.tau_maxvis)
        if chistar is not None:
            self.chistar = chistar
            
        self.zs      = np.linspace(0,max_zs,int(2401*max_zs/12.))
        self.chis    = self.results.comoving_radial_distance(self.zs)

        self.dchis   = (self.chis[2:]-self.chis[:-2])/2.
        self.dzs     = (self.zs[2:]-self.zs[:-2])/2.
        self.chis    = (self.chis[1:-1])
        self.zs      = (self.zs[1:-1])
        self.a       = 1./(1.+self.zs)
        self.dz      = self.zs[1]-self.zs[0]

        gethz        = np.vectorize(self.results.h_of_z)
        self.hz      = gethz(self.zs) # convert everything in units of Mpc/h

        self.chish    = self.chis*self.h # Mpc/h now
        self.chistarh = self.chistar*self.h
        self.dchish   = self.dchis*self.h

        #-----------------------------------------------------------------------
        #log version
        self.chistarlog = (self.results.conformal_time(0)- self.results.tau_maxvis)
        if chistar is not None:
            self.chistarlog = chistar
        self.chislog    = np.linspace(0,self.chistarlog,self.nz)
        self.zslog      = self.results.redshift_at_comoving_radial_distance(self.chislog)

        self.chislog    = self.results.comoving_radial_distance(self.zslog)
        self.dchislog   = (self.chislog[2:]-self.chislog[:-2])/2.
        self.dzslog     = (self.zslog[2:]-self.zslog[:-2])/2.
        self.chislog    = (self.chislog[1:-1])
        self.zslog      = (self.zslog[1:-1])
        self.alog       = 1./(1.+self.zslog)
        self.dzlog      = self.zslog[1]-self.zslog[0]

        gethzlog        = np.vectorize(self.results.h_of_z)
        self.hzlog      = gethz(self.zslog) # convert everything in units of Mpc/h

        self.chishlog    = self.chislog*self.h # Mpc/h now
        self.chistarhlog = self.chistarlog*self.h
        self.dchishlog   = self.dchislog*self.h

        if (feedback>2):
            table=[['get_lenspotential_cls()'],['get_lensed_scalar_cls()'],['get_unlensed_scalar_cls()'],[' ']] 
            print(tabulate(table,headers=['PrimaryCls']))   
            
            table=[['get_Wcmb()'],['get_Wcmblog()'],['get_Wgal(nz_lens)'],['get_Wshear(nz_src)'],['get_WIA()'],[' ']] 
            print(tabulate(table,headers=['Kernels'])) 

            table=[['get_nlnn(nbar_lens)'],['get_nlgg(nbar_src,sige)'],[' ']] 
            print(tabulate(table,headers=['Noise spectrum'])) 

    def get_lenspotential_cls(self,raw_cl=False,verbose=True,CMB_unit='muK',return_cls=False):
        if verbose==True:
            print("PP, PT, PE starting with ell=0")
        lens_potential_cls = self.results.get_lens_potential_cls(lmax=self.lmax,raw_cl=raw_cl,CMB_unit='muK')
        self.lens_potential_cls = lens_potential_cls
        if return_cls==True:
            return lens_potential_cls[:,0],lens_potential_cls[:,1],lens_potential_cls[:,2]

    def get_lensed_scalar_cls(self,raw_cl=False,verbose=True,CMB_unit='muK',return_cls=False):
        if verbose==True:
            print("TT, EE, BB, TE")
        lensed_scalar_cls  = self.results.get_lensed_scalar_cls(lmax=self.lmax,raw_cl=raw_cl,CMB_unit='muK')
        self.lensed_scalar_cls = lensed_scalar_cls
        if return_cls==True:
            return lensed_scalar_cls[:,0],lensed_scalar_cls[:,1],lensed_scalar_cls[:,2],lensed_scalar_cls[:,3]

    def get_unlensed_scalar_cls(self,raw_cl=False,verbose=True,CMB_unit='muK',return_cls=False):
        if verbose==True:
            print("TT, EE, BB, TE")
        unlensed_scalar_cls = self.results.get_unlensed_scalar_cls(lmax=self.lmax,raw_cl=raw_cl,CMB_unit='muK')
        self.unlensed_scalar_cls = unlensed_scalar_cls
        if return_cls==True:
            return unlensed_scalar_cls[:,0],unlensed_scalar_cls[:,1],unlensed_scalar_cls[:,2],unlensed_scalar_cls[:,3]

    def get_Wgal(self,nz_len,bias):
        
        Wgal = np.zeros( (len(self.zs),nz_len.shape[1]-1))
        
        for i in range(0,nz_len.shape[1]-1):
            nz        = np.interp(self.zs,nz_len[:,0],nz_len[:,i+1],left=0,right=0)
            Wgal[:,i] = nz*bias[i]*self.dzs/self.dchish
        self.Wgal  = Wgal


    def get_Wshear(self,nz_src):
        
        Wshear=np.zeros((self.zs.shape[0],nz_src.shape[1]-1))

        for j in range(0,nz_src.shape[1]-1):
            nz        = np.interp(self.zs,nz_src[:,0],nz_src[:,j+1],left=0,right=0)
            for i in range(0,len(self.zs)):
                idx, = np.where(self.zs>=self.zs[i])
                tmp  = np.sum(nz[idx]*self.dzs[idx]/self.dchish[idx]*(self.chish[idx] - self.chish[i])/self.chish[idx]*self.dchish[idx])
                Wshear[i,j] = self.fc*self.chish[i]/self.a[i]*tmp

        self.Wshear=Wshear

    def get_Wia_nla(self,nz_src,A_IA=0.44,eta_IA=-0.7,z0=0.62):
        
        nz  = np.interp(self.zs,nz_src[:,0],nz_src[:,j+1],left=0,right=0)

        kh, zt, pk = self.results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
        s8  = np.array(self.results.get_sigma8())
        s8f = np.flip(s8)
        Dz  = s8f/np.max(s8f)

        C   = 3*(100*1000/3.086e22)**2/6.67e-11/8/np.pi/(3.24e-23)**3/1.99e30*5e-14 #C1*rhocrit
        Wia = -A_IA*C*self.omm/np.interp(self.zs,zt,Dz,left=1e30,right=1e30)*((1+self.zs)/(1+z0))**eta_IA*nz*self.dzs/self.dchish

        self.Wia=Wia

    def get_Wcmb(self):
        print("appending Wcmb: Note this only goes up to z=12")
        print("For full CMB lensing kernel use get_Wcmblog()")
        Wcmb       = ((self.chistarh-self.chish)/(self.chish**2*self.chistarh))
        Wcmba      = self.fc/self.a*(self.chistarh-self.chish)/(self.chistarh)*self.chish
        self.Wcmb  = Wcmb
        self.Wcmba = Wcmba

    def get_Wcmblog(self):
        print("appending Wcmblog")
        Wcmblog       = ((self.chistarhlog-self.chishlog)/(self.chishlog**2*self.chistarhlog))
        Wcmbalog      = self.fc/self.alog*(self.chistarhlog-self.chishlog)/(self.chistarhlog)*self.chishlog
        self.Wcmblog  = Wcmblog  # Using Weyl 
        self.Wcmbalog = Wcmbalog # Not using Weyl

    def limber(self,xtype,zmax=0,nonlinear=True):
        
        cl = np.zeros(self.lmax+1)

        if xtype=='nn':
            transftype1 = model.Transfer_tot
            transftype2 = model.Transfer_tot
            nbinN1 = self.Wgal.shape[1]
            nbinN2 = self.Wgal.shape[1]
            print("nbins n1 = %d"%nbinN1)
            print("nbins n2 = %d"%nbinN2)
            tmp       = np.zeros((self.lmax+1)) 
            self.clnn = np.zeros((nbinN1,nbinN2,self.lmax+1)) 
            PK = camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=4.0,zmin=0.0,var1=transftype1,var2=transftype2)
            
            for n1 in range(0,nbinN1):
                for n2 in range(0,nbinN2):
                    for l in range(1,self.lmax+1):
                        k      = (l+0.5)/self.chish[:]
                        pp     = PK.P(self.zs,k,grid=False)
                        tmp[l] = np.dot(self.dchish[:]/self.chish[:]**2*self.Wgal[:,n1]*self.Wgal[:,n2],pp)
                    self.clnn[n1,n2,:]=tmp
                       

        if xtype=='gg':
            transftype1 = model.Transfer_tot
            transftype2 = model.Transfer_tot
            nbinG1      = self.Wshear.shape[1]
            nbinG2      = self.Wshear.shape[1]
            print("nbins g1 = %d"%nbinG1)
            print("nbins g2 = %d"%nbinG2)
            tmp  = np.zeros((self.lmax+1)) 
            self.clgg = np.zeros((nbinG1,nbinG2,self.lmax+1))  
            PK = camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=4.0,zmin=0.0,var1=transftype1,var2=transftype2)
            c=0
            for n1 in range(0,nbinG1):
                for n2 in range(0,nbinG2):
                    for l in frogress.bar(range(1,self.lmax+1)):
                        k  = (l+0.5)/self.chish[:]
                        pp  = PK.P(self.zs,k,grid=False)
                        tmp[l]=np.dot(self.dchish[:]/self.chish[:]**2*self.Wshear[:,n1]*self.Wshear[:,n2],pp)  
                    self.clgg[n1,n2,:]=tmp
                
        if xtype=='ng':
            transftype1 = model.Transfer_tot
            transftype2 = model.Transfer_tot
            nbin1       = self.Wgal.shape[1]
            nbin2       = self.Wshear.shape[1]
            print("nbins n1 = %d"%nbin1)
            print("nbins g2 = %d"%nbin2)
            tmp  = np.zeros((self.lmax+1)) 
            self.clng = np.zeros((nbin1,nbin2,self.lmax+1)) 
            PK = camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=4.0,zmin=0.0,var1=transftype1,var2=transftype2)
            for n1 in range(0,nbin1):
                for n2 in range(0,nbin2):
                    for l in range(1,self.lmax+1):
                        k      = (l+0.5)/self.chish[:]
                        pp     = PK.P(self.zs,k,grid=False)
                        tmp[l] = np.dot(self.dchish[:]/self.chish[:]**2*self.Wgal[:,n1]*self.Wshear[:,n2],pp)
                    self.clng[n1,n2,:]=tmp

        if xtype=='nk':
            transftype1 = model.Transfer_tot
            transftype2 = model.Transfer_tot
            nbin1       = self.Wgal.shape[1]
            nbin2       = 1
            print("nbins n1 = %d"%nbin1)
            print("nbins k1 = %d"%nbin2)
            tmp  = np.zeros((self.lmax+1)) 
            self.clnk = np.zeros((nbin1,nbin2,self.lmax+1)) 
            PK = camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=4.0,zmin=0.0,var1=transftype1,var2=transftype2)
            for n1 in range(0,nbin1):
                for l in range(1,self.lmax+1):
                    k  = (l+0.5)/self.chish[:]
                    pp  = PK.P(self.zs,k,grid=False)
                    tmp[l]=np.dot(self.dchish[:]/self.chish[:]**2*self.Wgal[:,n1]*self.Wcmba[:],pp)
                self.clnk[n1,0,:]=tmp

        if xtype=='gk':
            transftype1 = model.Transfer_tot
            transftype2 = model.Transfer_tot
            nbin1       = self.Wshear.shape[1]
            nbin2       = 1
            print("nbins g1 = %d"%nbin1)
            print("nbins k1 = %d"%nbin2)  
            tmp  = np.zeros((self.lmax+1)) 
            self.clgk = np.zeros((nbinG,1,self.lmax+1))  
            PK    = camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=4.0,zmin=0.0,var1=transftype1,var2=transftype2)
            for n1 in range(0,nbin1):
                for l in range(1,self.lmax+1):
                    k     = (l+0.5)/self.chish[:]
                    pp    = PK.P(self.zs,k,grid=False)
                    tmp[l] = np.dot(self.dchish[:]/self.chish[:]**2*self.Wshear[:,n1]*self.Wcmba,pp)
                self.clgk[n1,0,:] = tmp
                      
        if xtype=='kklog':
            transftype1 = model.Transfer_Weyl
            transftype2 = model.Transfer_Weyl
            print("nbins k1 = %d"%1)
            print("nbins k2 = %d"%1)
            self.clkk = np.zeros((1,1,self.lmax+1))

            PK  =  camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=1100,zmin=0.0,var1=transftype1,var2=transftype2)
            
            if zmax==0: zmax=1200
            idx = np.where( (self.zslog<zmax) )[0]
            
            for l in range(1,self.lmax+1):
                k     = (l+0.5)/self.chishlog[idx]
                pp    = PK.P(self.zslog[idx],k,grid=False)
                cl[l] = np.dot(self.dchishlog[idx]*self.Wcmblog[idx]*self.Wcmblog[idx]/k**4/self.h**4,pp)*(l*(l+1))**2

            self.clkk[0,0,:] = cl

        if xtype=='kklin':
            if zmax>12:
                sys.exit('zmax set to less than 12. Must use kklog for that')

            transftype1 = model.Transfer_Weyl
            transftype2 = model.Transfer_Weyl
            self.clkk = np.zeros((1,1,self.lmax+1))

            PK  =  camb.get_matter_power_interpolator(self.pars,nonlinear=nonlinear,hubble_units=True, k_hunit=True, kmax=500, zmax=1100,zmin=0.0,var1=transftype1,var2=transftype2)
            
            if zmax==0: zmax=1200
            else:
                print("computing up to z=%.2f"%zmax)

            idx = np.where( (self.zs<zmax) )[0]   
            
            for l in range(1,self.lmax+1):
                k     = (l+0.5)/self.chish[idx]
                pp    = PK.P(self.zs[idx],k,grid=False)
                cl[l] = np.dot(self.dchish[idx]*self.Wcmb[idx]*self.Wcmb[idx]/k**4/self.h**4,pp)*(l*(l+1))**2
                
                #cl[l] = np.dot(self.dchish[idx]/self.chish[idx]**2*self.Wcmba[idx]*self.Wcmba[idx],pp) * (l*(l+1))**2
                
                   

            self.clkk[0,0,:] = cl

