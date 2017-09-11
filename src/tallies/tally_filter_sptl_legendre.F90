module tally_filter_sptl_legendre

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants
  use error
  use hdf5_interface
  use math,                only: calc_pn
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! SptlLegendreFilter gives Legendre moments of the change in scattering angle
!===============================================================================

  type, public, extends(TallyFilter) :: SpatialLegendreFilter
    integer :: order
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type SpatialLegendreFilter

contains

!===============================================================================
! SpatialLegendreFilter methods
!===============================================================================

  subroutine from_xml(this, node)
    class(SpatialLegendreFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    ! Get specified order
    call get_node_value(node, "order", this % order)
    this % n_bins = this % order + 1
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(SpatialLegendreFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: i, j, n
    integer :: num_nm
    real(8) :: rn(2*this % order + 1)

    ! TODO: Use recursive formula to calculate higher orders
    j = 0
    do n = 0, this % order
      ! Calculate n-th order
      num_nm = 2*n + 1
      rn(1:num_nm) = calc_pn(n, p % last_xyz)

      ! Append matching (bin,weight) for each moment
      do i = 1, num_nm
        j = j + 1
        call match % bins % push_back(j)
        call match % weights % push_back(rn(i))
      end do
    end do
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(SpatialLegendreFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "spatiallegendre")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(SpatialLegendreFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    label = "Legendre expansion, P" // trim(to_str(bin - 1))
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================


end module tally_filter_sptl_legendre
