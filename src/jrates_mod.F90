!-------------------------------------------------------------------------------
! photolysis rate constants
!-------------------------------------------------------------------------------
module jrates_mod
  use phot_kind_mod, only: rk => kind_phot
  use dissociation_class, only: dissociation_t
  use o3_dissociation_class, only: o3_dissociation_t
  use phot_specs_mod, only: phot_specs_t, phot_specs_array, NPhotSpecifiers, phot_specs_init

  implicit none

  private
  public :: jrates_init
  public :: jrates_calc

  ! linked list of dissociation objects
  type dissListElem
     type(phot_specs_t), pointer :: rxn_specifier => null()
     class(dissociation_t), pointer :: diss_obj => null()
     type(dissListElem), pointer :: next => null()
  end type dissListElem
  type(dissListElem), pointer :: diss_list => null()

  ! used to map the diss objs back to the list of jnames
  type ratesMap
     integer :: jndx
     class(dissociation_t), pointer :: diss_obj => null()
     character(len=16) :: jname
     integer :: channel_number = -1
  end type ratesMap
  type(ratesMap), allocatable :: ratesMapArray(:)

contains

  !-----------------------------------------------------------------------------
  ! initialize photo-dissociation objects based on reaction tag names
  !-----------------------------------------------------------------------------
  subroutine  jrates_init( jnames, errmsg, errflg )
    use params_mod, only: input_data_root

    character(len=*), intent(in)  :: jnames(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    type(dissListElem), pointer :: ptr => null()
    integer :: i, nj, j, n

    ! the molecule dissociation objects (instantiated only if needed)
    type(o3_dissociation_t), pointer :: o3_diss_obj => null()
    type(dissociation_t), pointer :: diss_obj => null()

    call phot_specs_init(input_data_root)

    nj = size(jnames)
    allocate(ratesMapArray(nj))

    ! build a collection of dissociation objects
    do i = 1,NPhotSpecifiers
       findj: do j = 1,nj
          if ( any( phot_specs_array(i)%tag_names(:)==jnames(j) ) ) then
             select case (phot_specs_array(i)%diss_type)
             case ('dissociation_t')
                allocate(diss_obj)
                call add_to_list( phot_specs_array(i),diss_obj )
             case ('o3_dissociation_t')
                allocate(o3_diss_obj)
                call add_to_list( phot_specs_array(i),o3_diss_obj )
             case default
                errmsg = 'jrates_init: '//trim(jnames(j))//' not recognized'
                errflg = -1
                return
             end select
             exit findj
          end if
       end do findj
    end do

    ! initialize the dissociation objects
    ptr => diss_list
    do while (associated(ptr))
       call ptr%diss_obj%initialize(ptr%rxn_specifier, errmsg, errflg)
       ptr => ptr%next
    end do

    ! set up mapping of instantiated diss objs to desired Js
    do j = 1,nj
       ptr => diss_list
       do while (associated(ptr))
          if ( any( ptr%diss_obj%tag_names(:)==jnames(j) ) ) then
             ratesMapArray(j)%diss_obj=>ptr%diss_obj
             ratesMapArray(j)%jname = jnames(j)
             do n = 1,ratesMapArray(j)%diss_obj%nchannels
                if (ratesMapArray(j)%jname == ratesMapArray(j)%diss_obj%tag_names(n)) then
                   ratesMapArray(j)%channel_number = n
                end if
             end do
          end if
          ptr => ptr%next
       end do
       if ( .not. associated(ratesMapArray(j)%diss_obj)) then
          errflg = -1
          errmsg = 'jrates_init: not able to map '//trim(jnames(i))//' with a diss obj'
          return
       end if
    end do

  end subroutine jrates_init

  !-----------------------------------------------------------------------------
  ! calculate photolysis rate constants
  !-----------------------------------------------------------------------------
  subroutine jrates_calc( radfld, etf, esfact, temp, dens, rate_consts, rxn_eqns )

    real(rk), intent(in) :: radfld(:,:)
    real(rk), intent(in) :: etf(:)
    real(rk), intent(in) :: esfact
    real(rk), intent(in) :: temp(:)
    real(rk), intent(in) :: dens(:)

    real(rk), intent(out) :: rate_consts(:,:)
    character(len=*), optional, intent(out) :: rxn_eqns(:)

    integer :: j, nj, n
    type(dissListElem), pointer :: ptr

    nj = size(rate_consts,dim=2)

    ! update the photo-dissociation rate constants
    ptr => diss_list
    do while (associated(ptr))
       call ptr%diss_obj%update_rates(temp, dens, radfld, etf, esfact)
       ptr => ptr%next
    end do

    ! copy in the corresponding dissociation channel rate constants
    do j = 1,nj
       n = ratesMapArray(j)%channel_number
       rate_consts(:,j) = ratesMapArray(j)%diss_obj%rateconsts(:,n)
    end do

    ! copy out the reaction descriptors if requested
    if (present(rxn_eqns)) then
       do j = 1,nj
          n = ratesMapArray(j)%channel_number
          rxn_eqns(j) = ratesMapArray(j)%diss_obj%equations(n)
       end do
    end if

  end subroutine jrates_calc

  ! private routines
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  subroutine add_to_list( rxn_specifier, diss_obj )
    type(phot_specs_t), target,     intent(in) :: rxn_specifier
    class(dissociation_t), pointer, intent(in) :: diss_obj

    type(dissListElem), pointer :: newelem

    allocate(newelem)
    newelem%rxn_specifier => rxn_specifier
    newelem%diss_obj => diss_obj

    diss_list => insertList( diss_list, newelem )

  end subroutine add_to_list

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  function insertList( head, elem)

    type(dissListElem), pointer :: head, elem
    type(dissListElem), pointer :: insertList

    elem%next => head
    insertList=> elem

  end function insertList

end module jrates_mod
