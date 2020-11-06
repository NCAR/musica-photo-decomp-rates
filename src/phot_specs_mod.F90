module phot_specs_mod
  use phot_kind_mod, only: rk => kind_phot

  integer, parameter :: NPhotSpecifiers = 12

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

    integer :: m

    ! could read a configuration file
    ! for now just hard wire the specifiers

    m = 0
    allocate(phot_specs_array(NPhotSpecifiers))

    m=m+1
    phot_specs_array(m)%diss_type = 'o3_dissociation_t'
    phot_specs_array(m)%nchannels = 2
    phot_specs_array(m)%tag_names(1) = 'jo3a'
    phot_specs_array(m)%equations(1) = 'O3 + hv -> O2 + O(1D)'
    phot_specs_array(m)%tag_names(2) = 'jo3b'
    phot_specs_array(m)%equations(2) = 'O3 + hv -> O2 + O(3P)'
    phot_specs_array(m)%quantum_yields = -huge(1._rk)
    phot_specs_array(m)%nfiles = 3
    phot_specs_array(m)%nskip(1) = 9
    phot_specs_array(m)%nread(1) = 56
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_O3_Ackerman_1971.txt'
    phot_specs_array(m)%nskip(2) = 28
    phot_specs_array(m)%nread(2) = 65
    phot_specs_array(m)%filename(2) = trim(input_data_root)//'/XSQY/XS_O3_218_JPL06.txt'
    phot_specs_array(m)%nskip(3) = 39
    phot_specs_array(m)%nread(3) = 168
    phot_specs_array(m)%filename(3) = trim(input_data_root)//'/XSQY/XS_O3_298_JPL06.txt'

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 2
    phot_specs_array(m)%tag_names(1) = 'jch4a'
    phot_specs_array(m)%equations(1) = 'CH4 + hv -> CH3O2 + H'
    phot_specs_array(m)%tag_names(2) = 'jch4b'
    phot_specs_array(m)%equations(2) = 'CH4 + hv -> H2 + Products'
    phot_specs_array(m)%quantum_yields(:2) = (/0.45_rk,0.55_rk/)
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_CH4.txt'
    phot_specs_array(m)%nskip(1) = 22
    phot_specs_array(m)%nread(1) = 39
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'CO2 + hv -> CO + O'
    phot_specs_array(m)%tag_names(1) = 'jco2'
    phot_specs_array(m)%quantum_yields(:1) = (/1._rk/)
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_CO2.txt'
    phot_specs_array(m)%nskip(1) = 20
    phot_specs_array(m)%nread(1) = 55
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'COF2 + hv  -> 2F'
    phot_specs_array(m)%tag_names(1) = 'jcof2'
    phot_specs_array(m)%quantum_yields(:1) = (/1._rk/)
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_COF2_JPL10.txt'
    phot_specs_array(m)%nskip(1) = 31
    phot_specs_array(m)%nread(1) = 21
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'COFCl + hv  -> Cl + F'
    phot_specs_array(m)%tag_names(1) = 'jcofcl'
    phot_specs_array(m)%quantum_yields(:1) = (/1._rk/)
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_COFCL_JPL10.txt'
    phot_specs_array(m)%nskip(1) = 32
    phot_specs_array(m)%nread(1) = 32
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'HBr + hv -> H + Br'
    phot_specs_array(m)%tag_names(1) = 'jhbr'
    phot_specs_array(m)%quantum_yields(:1) = (/1._rk/)
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_HBR_JPL06.txt'
    phot_specs_array(m)%nskip(1) = 44
    phot_specs_array(m)%nread(1) = 40
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'HF + hv  -> H + F'
    phot_specs_array(m)%tag_names(1) = 'jhf'
    phot_specs_array(m)%quantum_yields(1) = 1._rk
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_HF.txt'
    phot_specs_array(m)%nskip(1) = 14
    phot_specs_array(m)%nread(1) = 39
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'SF6 + hv -> product'
    phot_specs_array(m)%tag_names(1) = 'jsf6'
    phot_specs_array(m)%quantum_yields(1) = 1._rk
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_SF6.txt'
    phot_specs_array(m)%nskip(1) = 14
    phot_specs_array(m)%nread(1) = 14
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'GLYALD + hv -> 2HO2 + CO + CH2O' !  GLYALD (HOCH2CHO) + hv -> 2*HO2 + CO + CH2O
    phot_specs_array(m)%tag_names(1) = 'jglyald'
    phot_specs_array(m)%quantum_yields(1) = 0.5_rk
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_GLYALD.txt'
    phot_specs_array(m)%nskip(1) = 15
    phot_specs_array(m)%nread(1) = 131
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'HYAC + hv -> CH3CO3 + HO2 + CH2O'
    phot_specs_array(m)%tag_names(1) = 'jhyac'
    phot_specs_array(m)%quantum_yields(1) = 0.65_rk
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_HYAC.txt'
    phot_specs_array(m)%nskip(1) = 8
    phot_specs_array(m)%nread(1) = 101
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'dissociation_t'
    phot_specs_array(m)%nchannels = 2
    phot_specs_array(m)%equations(1) = 'MACR + hv -> 1.34HO2 + 0.66MCO3 + 1.34CH2O+ CH3CO3'
    phot_specs_array(m)%tag_names(1) = 'jmacr_a'
    phot_specs_array(m)%quantum_yields(1) = 0.005_rk
    phot_specs_array(m)%equations(2) = 'MACR + hv -> 0.66OH + 1.34CO'
    phot_specs_array(m)%tag_names(2) = 'jmacr_b'
    phot_specs_array(m)%quantum_yields(2) = 0.005_rk
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/XSQY/XS_MACR_JPL06.txt'
    phot_specs_array(m)%nskip(1) = 20
    phot_specs_array(m)%nread(1) = 146
    phot_specs_array(m)%xfac(1)  = 1._rk

    m=m+1
    phot_specs_array(m)%diss_type = 'no2_dissociation_t'
    phot_specs_array(m)%nchannels = 1
    phot_specs_array(m)%equations(1) = 'NO2 -> NO + O(3P)'
    phot_specs_array(m)%tag_names(1) = 'jno2'
    phot_specs_array(m)%nfiles = 1
    phot_specs_array(m)%filename(1) = trim(input_data_root)//'/DATAJ1/YLD/NO2_jpl11.yld'
    phot_specs_array(m)%nskip(1) = 2
    phot_specs_array(m)%nread(1) = 25
    phot_specs_array(m)%xfac(1)  = 1._rk

    if (m /= NPhotSpecifiers) then
       write(*,*) 'ERROR phot_specs_init: m /= NPhotSpecifiers.. m = ',m, ' NPhotSpecifiers = ',NPhotSpecifiers
       stop 99
    end if

  end subroutine phot_specs_init

end module phot_specs_mod
