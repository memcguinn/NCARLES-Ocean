      function stokes_ker(sigma)
c
c ----------- evaluate kernel of the stokes integral
c             for the donelan spectral shape
c             the latter is converted into "sigma" space
c             careful with 2*pi factors
c
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats

      wave_spec =  (ann*grav_w*grav_w/(sigma_p*sigma**4))*
     +             exp(-bnn*(sigma_p/sigma)**4)
      stokes_ker = 2.0*(wave_spec*sigma**3)*
     +             exp(2.0*sigma*sigma*z_pt/grav_w)/grav_w
c
      return
      end
