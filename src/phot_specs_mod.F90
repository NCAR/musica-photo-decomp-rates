module phot_specs_mod
  use phot_kind_mod, only: rk => kind_phot

  integer, parameter :: NPhotSpecifiers=2

  integer, parameter :: max_files = 5
  integer, parameter :: max_rxns = 4

  type phot_specs_t
     character(len=80)  :: diss_type
     integer            :: nchannels
     character(len=80)  :: equations(max_rxns)
     character(len=16)  :: tag_names(max_rxns)
     real(rk)           :: quantum_yields(max_rxns)
     integer            :: nfiles
     integer            :: nskip(max_files)
     integer            :: nread(max_files)
     real(rk)           :: xfac(max_files)
     character(len=400) :: filename(max_files)
  end type phot_specs_t

  type(phot_specs_t), protected, allocatable :: phot_specs_array(:)

contains

  subroutine phot_specs_init( input_data_root )
    character(len=*), intent(in) :: input_data_root

    ! could read a configuration file
    ! for now just hard wire the specifiers

    allocate(phot_specs_array(NPhotSpecifiers))

    phot_specs_array(1)%diss_type = 'o3_dissociation_t'
    phot_specs_array(1)%nchannels = 2
    phot_specs_array(1)%tag_names(1) = 'jo3a'
    phot_specs_array(1)%equations(1) = 'O3 + hv -> O2 + O(1D)'
    phot_specs_array(1)%tag_names(2) = 'jo3b'
    phot_specs_array(1)%equations(2) = 'O3 + hv -> O2 + O(3P)'
    phot_specs_array(1)%quantum_yields = -huge(1._rk)
    phot_specs_array(1)%nfiles = 3
    phot_specs_array(1)%nskip(1) = 9
    phot_specs_array(1)%nread(1) = 56
    phot_specs_array(1)%filename(1) = trim(input_data_root)//'/XSQY/XS_O3_Ackerman_1971.txt'
    phot_specs_array(1)%nskip(2) = 28
    phot_specs_array(1)%nread(2) = 65
    phot_specs_array(1)%filename(2) = trim(input_data_root)//'/XSQY/XS_O3_218_JPL06.txt'
    phot_specs_array(1)%nskip(3) = 39
    phot_specs_array(1)%nread(3) = 168
    phot_specs_array(1)%filename(3) = trim(input_data_root)//'/XSQY/XS_O3_298_JPL06.txt'

    phot_specs_array(2)%diss_type = 'dissociation_t'
    phot_specs_array(2)%nchannels = 2
    phot_specs_array(2)%tag_names(1) = 'jch4a'
    phot_specs_array(2)%equations(1) = 'CH4 + hv -> CH3O2 + H'
    phot_specs_array(2)%tag_names(2) = 'jch4b'
    phot_specs_array(2)%equations(2) = 'CH4 + hv -> H2 + Products'
    phot_specs_array(2)%quantum_yields(:2) = (/0.45_rk,0.55_rk/)
    phot_specs_array(2)%nfiles = 1
    phot_specs_array(2)%filename(1) = trim(input_data_root)//'/XSQY/XS_CH4.txt'
    phot_specs_array(2)%nskip(1) = 22
    phot_specs_array(2)%nread(1) = 39
    phot_specs_array(2)%xfac(1)  = 1._rk


  end subroutine phot_specs_init

end module phot_specs_mod
