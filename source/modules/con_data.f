      module con_data
c ----------------------------------------------------------------------
        use pars, only : nscl
        type con_d
        sequence
        real ::  zo, vk, vkin, vk74, vk74in,
     +           grav, gcp, fcor, fcor_h, zi, pi2,
     +           batagk, t00, batag, vgcont, ugcont,
     +           cdbtm, dtdzf(nscl), dtjump, ugal, divgls,
     +           z1, utausv, xl, yl, zl, dx, dy, dz, dt,
     +           fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +           time, zody, zody74, tsfcc(nscl),
     +           utau, qstar(nscl), wtsfc(nscl),
     +           uwsfc, vwsfc, amonin,
     +           zol, hol, smal_e, sml_eg, rho_a, rho_w
        end type con_d
        type(con_d), target :: c_c
        real, pointer ::
     +           zo, vk, vkin, vk74, vk74in,
     +           grav, gcp, fcor, fcor_h, zi, pi2,
     +           batagk, t00, batag, vgcont, ugcont,
     +           cdbtm, dtdzf(:), dtjump, ugal, divgls,
     +           z1, utausv, xl, yl, zl, dx, dy, dz, dt,
     +           fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +           time, zody, zody74, tsfcc(:),
     +           utau, qstar(:), wtsfc(:),
     +           uwsfc, vwsfc, amonin,
     +           zol, hol, smal_e, sml_eg, rho_a, rho_w
c ----------------------------------------------------------------------
        type con_h
        sequence
        real ::  x_hurr, y_hurr, y_hurr_i, xpt_les, ypt_les,
     +           h_speed, t_move, u_10, v_10, cd_10, tau_x,
     +           tau_y, qw_tot_aw
        end type con_h
        type(con_h), target :: c_hurr
        real, pointer ::
     +           x_hurr, y_hurr, y_hurr_i, xpt_les, ypt_les,
     +           h_speed, t_move, u_10, v_10, cd_10, tau_x,
     +           tau_y, qw_tot_aw
      contains
         subroutine fill_cc
c
c --------------- pointer associations for constant variables
c
             zo     => c_c%zo
             vk     => c_c%vk
             vkin   => c_c%vkin
             vk74   => c_c%vk74
             vk74in => c_c%vk74in
             grav   => c_c%grav
             gcp    => c_c%gcp
             fcor   => c_c%fcor
             fcor_h => c_c%fcor_h
             zi     => c_c%zi
             pi2    => c_c%pi2
             batagk => c_c%batagk
             t00    => c_c%t00
             batag  => c_c%batag
             vgcont => c_c%vgcont
             ugcont => c_c%ugcont
             cdbtm  => c_c%cdbtm
             dtdzf  => c_c%dtdzf
             dtjump => c_c%dtjump
             ugal   => c_c%ugal
             divgls => c_c%divgls
             z1     => c_c%z1
             utausv => c_c%utausv
             xl     => c_c%xl
             yl     => c_c%yl
             zl     => c_c%zl
             dx     => c_c%dx
             dy     => c_c%dy
             dz     => c_c%dz
             dt     => c_c%dt
             fnxy   => c_c%fnxy
             dzdz   => c_c%dzdz
             dsl    => c_c%dsl
             c23    => c_c%c23
             dtgama => c_c%dtgama
             dtzeta => c_c%dtzeta
             xkmax  => c_c%xkmax
             time   => c_c%time
             zody   => c_c%zody
             zody74 => c_c%zody74
             tsfcc  => c_c%tsfcc
             utau   => c_c%utau
             qstar  => c_c%qstar
             wtsfc  => c_c%wtsfc
             uwsfc  => c_c%uwsfc
             vwsfc  => c_c%vwsfc
             amonin => c_c%amonin
             zol    => c_c%zol
             hol    => c_c%hol
             smal_e => c_c%smal_e
             sml_eg => c_c%sml_eg
             rho_a  => c_c%rho_a
             rho_w  => c_c%rho_w
         return
         end subroutine fill_cc
         subroutine fill_ch
c
c --------------- pointer associations for constant variables
c                 with hurricane information
c
             x_hurr    => c_hurr%x_hurr
             y_hurr    => c_hurr%y_hurr
             y_hurr_i  => c_hurr%y_hurr_i
             xpt_les   => c_hurr%xpt_les
             ypt_les   => c_hurr%ypt_les
             h_speed   => c_hurr%h_speed
             t_move    => c_hurr%t_move
             u_10      => c_hurr%u_10
             v_10      => c_hurr%v_10
             cd_10     => c_hurr%cd_10
             tau_x     => c_hurr%tau_x
             tau_y     => c_hurr%tau_y
             qw_tot_aw => c_hurr%qw_tot_aw
         return
         end subroutine fill_ch
      end module con_data
