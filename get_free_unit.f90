SUBROUTINE get_free_unit( funit )
  ! finds an available unit number for input/output
  ! avoids all units numbers less than 10, as some have special meaning.

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: funit
  INTEGER, SAVE :: last_funit = 10 ! avoid unit<10
  LOGICAL :: used
  
  funit = last_funit
  DO     
     INQUIRE( UNIT=funit, OPENED=used )
     IF ( .NOT. used ) EXIT
     funit = funit + 1
  END DO
  
  last_funit = funit
  
END SUBROUTINE get_free_unit
