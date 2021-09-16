!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qesgipaw_read_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes_1.0
  !
  USE FoX_dom
  USE qesgipaw_types_module
  !
  IMPLICIT NONE
  !
  PUBLIC qesgipaw_read
  !
  INTERFACE qesgipaw_read
    MODULE PROCEDURE qes_read_gipaw
    MODULE PROCEDURE qes_read_general_info
    MODULE PROCEDURE qes_read_format
    MODULE PROCEDURE qes_read_creator
    MODULE PROCEDURE qes_read_created
    MODULE PROCEDURE qes_read_parallel_info
    MODULE PROCEDURE qes_read_closed
    MODULE PROCEDURE qes_read_inputGIPAW
    MODULE PROCEDURE qes_read_job
    MODULE PROCEDURE qes_read_files
    MODULE PROCEDURE qes_read_scf_gipaw
    MODULE PROCEDURE qes_read_nmr
    MODULE PROCEDURE qes_read_efg
    MODULE PROCEDURE qes_read_hfi
  END INTERFACE qesgipaw_read
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_read_gipaw(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(gipaw_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "general_info")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gipawType","general_info: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gipawType","general_info: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%general_info_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_general_info(tmp_node, obj%general_info, ierr )
    ELSE
       obj%general_info_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "parallel_info")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gipawType","parallel_info: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gipawType","parallel_info: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%parallel_info_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_parallel_info(tmp_node, obj%parallel_info, ierr )
    ELSE
       obj%parallel_info_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "inputGIPAW")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gipawType","inputGIPAW: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gipawType","inputGIPAW: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_inputGIPAW(tmp_node, obj%inputGIPAW, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "status")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gipawType","status: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gipawType","status: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%status_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%status , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gipawType","error reading status")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gipawType","error reading status",10)
         END IF
      END IF
    ELSE
       obj%status_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "cputime")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gipawType","cputime: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gipawType","cputime: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%cputime_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%cputime , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gipawType","error reading cputime")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gipawType","error reading cputime",10)
         END IF
      END IF
    ELSE
       obj%cputime_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "closed")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gipawType","closed: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gipawType","closed: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%closed_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_closed(tmp_node, obj%closed, ierr )
    ELSE
       obj%closed_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_gipaw
  !
  !
  SUBROUTINE qes_read_general_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(general_info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "format")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","format: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","format: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_format(tmp_node, obj%format, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "creator")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","creator: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","creator: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_creator(tmp_node, obj%creator, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "created")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","created: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","created: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_created(tmp_node, obj%created, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "job")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","job: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","job: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%job, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:general_infoType","error reading job")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:general_infoType","error reading job",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_general_info
  !
  !
  SUBROUTINE qes_read_format(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(format_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !

    IF (hasAttribute(xml_node, "NAME")) THEN
      CALL extractDataAttribute(xml_node, "NAME", obj%NAME)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: formatType",&
                        "required attribute NAME not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: formatType",&
                      "required attribute NAME not found", 10 )
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "VERSION")) THEN
      CALL extractDataAttribute(xml_node, "VERSION", obj%VERSION)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: formatType",&
                        "required attribute VERSION not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: formatType",&
                      "required attribute VERSION not found", 10 )
      END IF
    END IF
    !



    !
    !
    CALL extractDataContent(xml_node, obj%format )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_format
  !
  !
  SUBROUTINE qes_read_creator(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(creator_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !

    IF (hasAttribute(xml_node, "NAME")) THEN
      CALL extractDataAttribute(xml_node, "NAME", obj%NAME)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: creatorType",&
                        "required attribute NAME not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: creatorType",&
                      "required attribute NAME not found", 10 )
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "VERSION")) THEN
      CALL extractDataAttribute(xml_node, "VERSION", obj%VERSION)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: creatorType",&
                        "required attribute VERSION not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: creatorType",&
                      "required attribute VERSION not found", 10 )
      END IF
    END IF
    !



    !
    !
    CALL extractDataContent(xml_node, obj%creator )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_creator
  !
  !
  SUBROUTINE qes_read_created(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(created_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !

    IF (hasAttribute(xml_node, "DATE")) THEN
      CALL extractDataAttribute(xml_node, "DATE", obj%DATE)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: createdType",&
                        "required attribute DATE not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: createdType",&
                      "required attribute DATE not found", 10 )
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "TIME")) THEN
      CALL extractDataAttribute(xml_node, "TIME", obj%TIME)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: createdType",&
                        "required attribute TIME not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: createdType",&
                      "required attribute TIME not found", 10 )
      END IF
    END IF
    !



    !
    !
    CALL extractDataContent(xml_node, obj%created )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_created
  !
  !
  SUBROUTINE qes_read_parallel_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(parallel_info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "nprocs")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nprocs: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nprocs: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nprocs, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading nprocs")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading nprocs",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nthreads")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nthreads: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nthreads: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nthreads, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading nthreads")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading nthreads",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ntasks")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","ntasks: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","ntasks: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ntasks, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading ntasks")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading ntasks",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nbgrp")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nbgrp: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nbgrp: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nbgrp, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading nbgrp")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading nbgrp",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "npool")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","npool: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","npool: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%npool, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading npool")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading npool",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ndiag")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","ndiag: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","ndiag: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ndiag, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading ndiag")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading ndiag",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nimages")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nimages: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nimages: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nimages_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nimages , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:parallel_infoType","error reading nimages")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:parallel_infoType","error reading nimages",10)
         END IF
      END IF
    ELSE
       obj%nimages_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_parallel_info
  !
  !
  SUBROUTINE qes_read_closed(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(closed_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !

    IF (hasAttribute(xml_node, "DATE")) THEN
      CALL extractDataAttribute(xml_node, "DATE", obj%DATE)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: closedType",&
                        "required attribute DATE not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: closedType",&
                      "required attribute DATE not found", 10 )
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "TIME")) THEN
      CALL extractDataAttribute(xml_node, "TIME", obj%TIME)
    ELSE
      IF ( PRESENT(ierr) ) THEN
         CALL infomsg ( "qes_read: closedType",&
                        "required attribute TIME not found" )
         ierr = ierr + 1
      ELSE
         CALL errore ("qes_read: closedType",&
                      "required attribute TIME not found", 10 )
      END IF
    END IF
    !



    !
    !
    CALL extractDataContent(xml_node, obj%closed )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_closed
  !
  !
  SUBROUTINE qes_read_inputGIPAW(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(inputGIPAW_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "control_job")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputGIPAWType","control_job: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputGIPAWType","control_job: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%control_job_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_job(tmp_node, obj%control_job, ierr )
    ELSE
       obj%control_job_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "scf_gipaw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputGIPAWType","scf_gipaw: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputGIPAWType","scf_gipaw: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%scf_gipaw_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_scf_gipaw(tmp_node, obj%scf_gipaw, ierr )
    ELSE
       obj%scf_gipaw_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "files")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputGIPAWType","files: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputGIPAWType","files: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%files_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_files(tmp_node, obj%files, ierr )
    ELSE
       obj%files_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nmr_input")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputGIPAWType","nmr_input: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputGIPAWType","nmr_input: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nmr_input_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_nmr(tmp_node, obj%nmr_input, ierr )
    ELSE
       obj%nmr_input_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "efg_input")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputGIPAWType","efg_input: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputGIPAWType","efg_input: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%efg_input_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_efg(tmp_node, obj%efg_input, ierr )
    ELSE
       obj%efg_input_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "hfi_input")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputGIPAWType","hfi_input: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputGIPAWType","hfi_input: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%hfi_input_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_hfi(tmp_node, obj%hfi_input, ierr )
    ELSE
       obj%hfi_input_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_inputGIPAW
  !
  !
  SUBROUTINE qes_read_job(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(job_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "job")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:jobType","job: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:jobType","job: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%job_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%job , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:jobType","error reading job")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:jobType","error reading job",10)
         END IF
      END IF
    ELSE
       obj%job_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "restart_mode")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:jobType","restart_mode: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:jobType","restart_mode: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%restart_mode_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%restart_mode , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:jobType","error reading restart_mode")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:jobType","error reading restart_mode",10)
         END IF
      END IF
    ELSE
       obj%restart_mode_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "verbosity")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:jobType","verbosity: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:jobType","verbosity: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%verbosity_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%verbosity , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:jobType","error reading verbosity")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:jobType","error reading verbosity",10)
         END IF
      END IF
    ELSE
       obj%verbosity_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "max_seconds")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:jobType","max_seconds: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:jobType","max_seconds: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%max_seconds_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%max_seconds , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:jobType","error reading max_seconds")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:jobType","error reading max_seconds",10)
         END IF
      END IF
    ELSE
       obj%max_seconds_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_job
  !
  !
  SUBROUTINE qes_read_files(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(files_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "prefix")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:filesType","prefix: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:filesType","prefix: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%prefix_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%prefix , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:filesType","error reading prefix")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:filesType","error reading prefix",10)
         END IF
      END IF
    ELSE
       obj%prefix_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tmp_dir")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:filesType","tmp_dir: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:filesType","tmp_dir: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%tmp_dir_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%tmp_dir , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:filesType","error reading tmp_dir")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:filesType","error reading tmp_dir",10)
         END IF
      END IF
    ELSE
       obj%tmp_dir_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "filcurr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:filesType","filcurr: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:filesType","filcurr: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%filcurr_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%filcurr , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:filesType","error reading filcurr")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:filesType","error reading filcurr",10)
         END IF
      END IF
    ELSE
       obj%filcurr_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fildfield")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:filesType","fildfield: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:filesType","fildfield: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fildfield_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fildfield , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:filesType","error reading fildfield")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:filesType","error reading fildfield",10)
         END IF
      END IF
    ELSE
       obj%fildfield_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "filnics")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:filesType","filnics: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:filesType","filnics: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%filnics_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%filnics , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:filesType","error reading filnics")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:filesType","error reading filnics",10)
         END IF
      END IF
    ELSE
       obj%filnics_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_files
  !
  !
  SUBROUTINE qes_read_scf_gipaw(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(scf_gipaw_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "diagonalization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_gipawType","diagonalization: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_gipawType","diagonalization: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%diagonalization_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%diagonalization , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:scf_gipawType","error reading diagonalization")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:scf_gipawType","error reading diagonalization",10)
         END IF
      END IF
    ELSE
       obj%diagonalization_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "conv_threshold")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_gipawType","conv_threshold: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_gipawType","conv_threshold: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%conv_threshold_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%conv_threshold , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:scf_gipawType","error reading conv_threshold")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:scf_gipawType","error reading conv_threshold",10)
         END IF
      END IF
    ELSE
       obj%conv_threshold_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "q_gipaw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_gipawType","q_gipaw: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_gipawType","q_gipaw: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%q_gipaw_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%q_gipaw , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:scf_gipawType","error reading q_gipaw")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:scf_gipawType","error reading q_gipaw",10)
         END IF
      END IF
    ELSE
       obj%q_gipaw_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "r_rand")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_gipawType","r_rand: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_gipawType","r_rand: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%r_rand_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%r_rand , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:scf_gipawType","error reading r_rand")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:scf_gipawType","error reading r_rand",10)
         END IF
      END IF
    ELSE
       obj%r_rand_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spline_ps")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_gipawType","spline_ps: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_gipawType","spline_ps: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%spline_ps_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%spline_ps , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:scf_gipawType","error reading spline_ps")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:scf_gipawType","error reading spline_ps",10)
         END IF
      END IF
    ELSE
       obj%spline_ps_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "paw_proj")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_gipawType","paw_proj: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_gipawType","paw_proj: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%paw_proj_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%paw_proj , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:scf_gipawType","error reading paw_proj")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:scf_gipawType","error reading paw_proj",10)
         END IF
      END IF
    ELSE
       obj%paw_proj_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_scf_gipaw
  !
  !
  SUBROUTINE qes_read_nmr(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(nmr_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "use_nmr_macroscopic_shape")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:nmrType","use_nmr_macroscopic_shape: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:nmrType","use_nmr_macroscopic_shape: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%use_nmr_macroscopic_shape_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%use_nmr_macroscopic_shape , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:nmrType","error reading use_nmr_macroscopic_shape")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:nmrType","error reading use_nmr_macroscopic_shape",10)
         END IF
      END IF
    ELSE
       obj%use_nmr_macroscopic_shape_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nmr_macroscopic_shape")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:nmrType","nmr_macroscopic_shape: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:nmrType","nmr_macroscopic_shape: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nmr_macroscopic_shape_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nmr_macroscopic_shape , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:nmrType","error reading nmr_macroscopic_shape")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:nmrType","error reading nmr_macroscopic_shape",10)
         END IF
      END IF
    ELSE
       obj%nmr_macroscopic_shape_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_nmr
  !
  !
  SUBROUTINE qes_read_efg(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(efg_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "efg_q")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:efgType","efg_q: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:efgType","efg_q: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%efg_q_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%efg_q , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:efgType","error reading efg_q")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:efgType","error reading efg_q",10)
         END IF
      END IF
    ELSE
       obj%efg_q_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_efg
  !
  !
  SUBROUTINE qes_read_hfi(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(hfi_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(INOUT)                  :: ierr
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !




    !
    tmp_node_list => getElementsByTagname(xml_node, "hfi_output_unit")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hfiType","hfi_output_unit: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hfiType","hfi_output_unit: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%hfi_output_unit_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%hfi_output_unit , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:hfiType","error reading hfi_output_unit")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:hfiType","error reading hfi_output_unit",10)
         END IF
      END IF
    ELSE
       obj%hfi_output_unit_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "hfi_nuclear_g_tensor")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hfiType","hfi_nuclear_g_tensor: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hfiType","hfi_nuclear_g_tensor: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%hfi_nuclear_g_tensor_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%hfi_nuclear_g_tensor , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:hfiType","error reading hfi_nuclear_g_tensor")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:hfiType","error reading hfi_nuclear_g_tensor",10)
         END IF
      END IF
    ELSE
       obj%hfi_nuclear_g_tensor_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "core_relax_method")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hfiType","core_relax_method: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hfiType","core_relax_method: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%core_relax_method_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%core_relax_method , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:hfiType","error reading core_relax_method")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:hfiType","error reading core_relax_method",10)
         END IF
      END IF
    ELSE
       obj%core_relax_method_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "hfi_via_reconstruction_only")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hfiType","hfi_via_reconstruction_only: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hfiType","hfi_via_reconstruction_only: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%hfi_via_reconstruction_only_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%hfi_via_reconstruction_only , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:hfiType","error reading hfi_via_reconstruction_only")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:hfiType","error reading hfi_via_reconstruction_only",10)
         END IF
      END IF
    ELSE
       obj%hfi_via_reconstruction_only_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_hfi
  !
  !
END MODULE qesgipaw_read_module