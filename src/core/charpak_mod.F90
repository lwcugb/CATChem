!> \file charpak_mod.F90
!! \brief Module for character string manipulation
!! \ingroup core_modules
!!
!! \author Original GEOS-Chem Development Team
!! \author CATChem Development Team
!! \date 2023
!!
!! This module provides utilities for character string manipulation and text
!! processing operations used throughout the CATChem atmospheric chemistry model.
!!
!! Original code from https://github.com/geoschem/geos-chem/
!!
!! \details
!! The character package module provides comprehensive string manipulation
!! functions including case conversion, text cleaning, splitting, replacement,
!! and formatted output operations.
!!
MODULE Charpak_Mod
    !
    ! !USES:
    !
      IMPLICIT NONE
      PRIVATE
    !
    ! !PUBLIC MEMBER FUNCTIONS:
    !
      PUBLIC  :: CleanText
      PUBLIC  :: CntMat
      PUBLIC  :: CopyTxt
      PUBLIC  :: CStrip
      PUBLIC  :: IsDigit
      PUBLIC  :: ReadOneLine
      PUBLIC  :: StrRepl
      PUBLIC  :: StrSplit
      PUBLIC  :: StrSqueeze
      PUBLIC  :: To_UpperCase
      PUBLIC  :: TranLc
      PUBLIC  :: TranUc
      PUBLIC  :: Txtext
      PUBLIC  :: WordWrapPrint
      PUBLIC  :: Unique
      PUBLIC  :: charArr2Str
      PUBLIC  :: str2CharArr
    !
    ! !REVISION HISTORY:
    !  See https://github.com/ufs-community/catchem                                  ! for complete history
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !DEFINED PARAMETERS
    !
      ! Maximum string length
      INTEGER, PARAMETER, PUBLIC :: MAXSTRLEN = 500

    CONTAINS

    ! \brief Count the number of times a character appears in a string
    !!
    !! \param  str1  String to search
    !! \param  str2  Character to search for
    !! \param imat  Number of matches
    !! \param locations  Positions of matches
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE CntMat( Str1, Str2, Imat, Locations )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN) ::  Str1             ! Text to scan
        CHARACTER(LEN=*), INTENT(IN) ::  Str2             ! Character to match
    !
    ! !OUTPUT PARAMETERS:
    !
        INTEGER,          INTENT(OUT) :: imat             ! Number of matches
        INTEGER,          OPTIONAL    :: Locations(255)   ! Positions of matches
    !
    ! !LOCAL VARIABLES:
    !
        ! Scalars
        INTEGER :: L1, L2, i, j
        LOGICAL :: again

        ! Arrays
        INTEGER :: TmpLocations(255)

        ! Initialize
        TmpLocations = 0
        L1           = MAX(1,LEN_TRIM(str1))
        L2           = LEN(str2)
        imat         = 0

        DO i=1,L1
           again = .true.
           j = 1
           DO WHILE (again)
              IF (str2(j:j).EQ.str1(i:i)) THEN
                 imat               = imat+1
                 TmpLocations(imat) = i
                 again              = .false.
              ELSEIF (j.LT.L2) THEN
                 j=j+1
              ELSE
                 again = .false.
              ENDIF
           ENDDO
        ENDDO

        ! Return positions where matches occurred (OPTIONAL)
        IF ( PRESENT( Locations ) ) Locations = TmpLocations

      END SUBROUTINE CntMat

      ! \brief Copy characters from one string to another
      !!
      !! \param  col  Starting column
      !! \param  str1  String to copy from
      !! \param str2  String to copy to
      !!
      !! \ingroup core_modules
      !!!>
      SUBROUTINE CopyTxt( col, str1, str2 )
    !
    ! !INPUT PARAMETERS:
    !
        INTEGER,          INTENT(IN)    :: col
        CHARACTER(LEN=*), INTENT(IN)    :: str1
    !
    ! !OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(INOUT) :: str2
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER :: ilt1,i1,i,j,ic

        i1 = LEN(str2)
        IF (i1.GT.0) THEN
           ilt1 = LEN(str1)
           IF (ilt1.GT.0) THEN
              ic = MAX0(col,1)
              i = 1
              j = ic
              DO WHILE ((i.LE.ilt1).and.(j.LE.i1))
                 str2(j:j) = str1(i:i)
                 i = i + 1
                 j = ic + (i-1)
              ENDDO
           ENDIF
        ENDIF

      END SUBROUTINE CopyTxt

    ! \brief Strip leading and trailing spaces from a string
    !!
    !! \param text        String to be modified
    !! \param    KeepSpaces  If =T, then keep spaces
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE CStrip( text, KeepSpaces )
    !
    ! !INPUT PARAMETERS:
    !
        LOGICAL,          OPTIONAL      :: KeepSpaces ! If =T, then keep spaces
                                                      !  but skip all other
                                                      !  non-printing chars
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(INOUT) :: TEXT       ! Text to be modified
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER          :: ilen, iasc, icnt, i, Start
        CHARACTER(LEN=1) :: ch

        ! Default: Skip space characters
        Start = 32

        ! If KEEPSPACES=T then skip all non-printing characters,
        ! but keep space characters. (bmy, 1/30/18)
        IF ( PRESENT( KeepSpaces ) ) THEN
           IF ( KeepSpaces ) Start = 31
        ENDIF

        ilen = LEN(text)
        IF (ilen.GT.1) THEN
           icnt = 1
           DO i=1,ilen
              iasc = ICHAR(text(i:i))

              ! Keep characters between these limits
              IF ( ( iasc > Start ).AND. (iasc < 255 ) ) THEN
                 ch = text(i:i)
                 text(icnt:icnt) = ch
                 icnt = icnt + 1
              ENDIF
           ENDDO
           ! Fill remainder of text with blanks
           DO i=icnt,ilen
              text(i:i) = ' '
           ENDDO
        ENDIF

      END SUBROUTINE CStrip

    ! \brief Check if a character is a digit
    !!
    !! \param  ch  Character to check
    !!
    !! \ingroup core_modules
    !!!>
      FUNCTION IsDigit( ch ) RESULT( lnum )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=1), INTENT(IN) :: ch
    !
    ! !RETURN VALUE:
    !
        LOGICAL                      :: lnum
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER iasc

        iasc = ICHAR(ch)
        lnum = .FALSE.
        IF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
           lnum = .TRUE.
        ENDIF

      END FUNCTION IsDigit
    ! \brief Replace text in a string
    !!
    !! \param str  String to be modified
    !! \param    Pattern  Pattern to search for
    !! \param    ReplTxt  Text to replace
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE StrRepl( Str, Pattern, ReplTxt )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN)    :: Pattern   ! Pattern to search for
        CHARACTER(LEN=*), INTENT(IN)    :: ReplTxt   ! Text to replace
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(INOUT) :: Str       ! String to be manipulated
    !
    ! !LOCAL VARIABLES:
    !
        ! Local variables
        INTEGER :: I1, I2

        !=================================================================
        ! StrRepl begins here!
        !=================================================================
        DO

           ! I1 is the first character that matches the search pattern;
           ! it must be 1 or larger.  Otherwise exit the routine.
           I1 = INDEX( Str, Pattern )
           IF ( I1 < 1 ) RETURN

           ! Replace the text.  I2 is the starting position of the
           ! string following the point of text replacement.
           I2 = I1 + LEN( Pattern )
           Str = Str(1:I1-1) // ReplTxt // Str(I2:)

        ENDDO

      END SUBROUTINE StrRepl
    ! \brief Split a string into substrings
    !!
    !! \param  Str  String to be searched
    !! \param  Sep  Separator character
    !! \param Result  Substrings
    !! \param N_SubStrs  # of substrings
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE StrSplit( Str, Sep, Result, N_SubStrs )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN)  :: STR           ! String to be searched
        CHARACTER(LEN=1), INTENT(IN)  :: SEP           ! Separator character
    !
    ! !OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(OUT) :: Result(255)   ! Returned substrings
        INTEGER,          OPTIONAL    :: N_SubStrs     ! # of substrings
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER             :: I, IFLAG, COL
        CHARACTER(LEN=2047) :: WORD

        !=======================================================================
        ! STRSPLIT begins here!
        !=======================================================================

        ! Initialize
        I         = 0
        COL       = 1
        IFLAG     = 0
        RESULT(:) = ''

        ! Loop until all matches found, or end of string
        DO WHILE ( IFLAG == 0 )

           ! Look for strings between separator string
           CALL TXTEXT ( SEP, TRIM( STR ), COL, WORD, IFLAG )

           ! Store substrings in RESULT array
           I         = I + 1
           RESULT(I) = TRIM( WORD )

        ENDDO

        ! Optional argument: return # of substrings found
        IF ( PRESENT( N_SUBSTRS ) ) N_SUBSTRS = I

      END SUBROUTINE StrSplit
    ! \brief Remove leading and trailing blanks from a string
    !!
    !! \param  Str  String to be manipulated
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE StrSqueeze( Str )
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(INOUT) :: Str   ! String to be squeezed

        !=================================================================
        ! STRSQUEEZE begins here!
        !=================================================================
        Str = ADJUSTR( TRIM( Str ) )
        Str = ADJUSTL( TRIM( Str ) )

      END SUBROUTINE StrSqueeze
    ! \brief Translate a character variable to all lower case letters.
    !!
    !! \param  text  String to be manipulated
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE TranLc( text )
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*) :: text
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER :: iasc,i,ilen

        ilen = LEN(text)
        DO I=1,ilen
           iasc = ICHAR(text(i:i))
           IF ((iasc.GT.64).AND.(iasc.LT.91)) THEN
              text(i:i) = CHAR(iasc+32)
           ENDIF
        ENDDO

      END SUBROUTINE TRANLC
    ! \brief Translate a character variable to all upper case letters.
    !!
    !! \param  text  String to be manipulated
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE TranUc( text )
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*) :: text
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER :: iasc,i,ilen

        ilen = LEN(text)
        DO i=1,ilen
           iasc = ICHAR(text(i:i))
           IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
              text(i:i) = CHAR(iasc-32)
           ENDIF
        ENDDO

      END SUBROUTINE TRANUC
    ! \brief Extract text from a string
    !!
    !! \param  ch  Character string
    !! \param  text  String to be manipulated
    !! \param  col  Column position
    !! \param  word  Extracted text
    !! \param  iflg  0: normal termination; 1: text not found
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE TxtExt(ch,text,col,word,iflg)
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN)    :: ch,text
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
        INTEGER,          INTENT(INOUT) :: col
    !
    ! !OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(OUT)   :: word
        INTEGER                         :: iflg
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER :: Tmax,T1,T2,imat
        LOGICAL :: again,prev

        ! Length of text
        Tmax = LEN(text)

        ! Fill Word with blanks
        WORD = REPEAT( ' ', LEN( WORD ) )

        IF (col.GT.Tmax) THEN
           ! Text does not contain any characters past Tmax.
           ! Reset col to one and return flag = {error condition}
           iflg = 1
           col = 1
        ELSEIF (col.EQ.Tmax) THEN
           ! End of TEXT reached
           CALL CntMat(ch,text(Tmax:Tmax),imat)
           IF (imat.EQ.0) THEN
              ! Copy character into Word and set col=1
              CALL CopyTxt(1,Text(Tmax:Tmax),Word)
              col = 1
              iflg = -1
           ELSE
              ! Same error condition as if col.GT.Tmax
              iflg = 1
           ENDIF
        ELSE
           ! Make sure column is not less than 1
           IF (col.LT.1) col=1
           CALL CntMat(ch,text(col:col),imat)
           IF (imat.GT.0) THEN
              prev=.true.
           ELSE
              prev=.false.
           ENDIF
           T1=col
           T2 = T1

           again = .true.
           DO WHILE (again)
              ! Check for a match with a character in ch
              CALL CntMat(ch,text(T2:T2),imat)
              IF (imat.GT.0) THEN
                 ! Current character in TEXT matches one (or more) of the
                 ! characters in ch.
                 IF (prev) THEN
                    IF (T2.LT.Tmax) THEN
                       ! Keep searching for a block of text
                       T2=T2+1
                       T1=T2
                    ELSE
                       ! Did not find any text blocks before running
                       ! out of characters in TEXT.
                       again=.false.
                       iflg=1
                    ENDIF
                 ELSE
                     ! Previous character did not match ch, so terminate.
                     ! NOTE: This is "NORMAL" termination of the loop
                    again=.false.
                    T2=T2-1
                    iflg = 0
                 ENDIF
              ELSEIF (T2.LT.Tmax) THEN
                 ! Add a letter to the current block of text
                 prev = .false.
                 T2=T2+1
              ELSE
                 ! Reached the end of the characters in TEXT before reaching
                 ! another delimiting character.  A text block was identified
                 ! however.
                 again=.false.
                 iflg=-1
              ENDIF
           ENDDO

           IF (iflg.EQ.0) THEN
              ! Copy characters into WORD and set col for return
              CALL CopyTxt(1,Text(T1:T2),Word)
              col = T2+1
           ELSE
              ! Copy characters into WORD and set col for return
              CALL CopyTxt(1,Text(T1:T2),Word)
              col = 1
           ENDIF
        ENDIF

      END SUBROUTINE TxtExt
    ! \brief Convert a string to uppercase
    !!
    !! \param  Text  String to be converted to uppercase
    !! \return  UpCaseText  Uppercase version of the input string
    !!
    !! \ingroup core_modules
    !!!>
      FUNCTION To_UpperCase( Text ) RESULT( UpCaseText )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN) :: Text         ! Input test
    !
    ! !RETURN VALUE:
    !
        CHARACTER(LEN=255)           :: UpCaseText   ! Output text, uppercase
    !
    ! !LOCAL VARIABLES:
    !
        ! Scalars
        INTEGER :: C, Ascii

        !=======================================================================
        ! Convert to uppercase
        !=======================================================================

        ! Initialize
        UpCaseText = Text

        ! Loop over all characters
        DO C = 1, LEN_TRIM( UpCaseText )

           ! Get the ASCII code for each character
           Ascii = ICHAR( UpCaseText(C:C) )

           ! If lowercase, convert to uppercase
           IF ( Ascii > 96 .and. Ascii < 123 ) THEN
              UpCaseText(C:C) = CHAR( Ascii - 32 )
           ENDIF
        ENDDO

      END FUNCTION To_UpperCase
    ! \brief Read one line from a file
    !!
    !! \param  fId    File unit number
    !! \param EndOfFile  Denotes EOF condition
    !! \param IoStatus   I/O status code
    !! \param  Squeeze    Call Strsqueeze?
    !!
    !! \ingroup core_modules
    !! \return Line read from file
    !!!>
      FUNCTION ReadOneLine( fId, EndOfFile, IoStatus, Squeeze ) RESULT( Line )
    !
    ! !INPUT PARAMETERS:
    !
        INTEGER, INTENT(IN)      :: fId        ! File unit number
        LOGICAL, OPTIONAL        :: Squeeze    ! Call Strsqueeze?
    !
    ! !OUTPUT PARAMETERS:
    !
        LOGICAL, INTENT(OUT)     :: EndOfFile  ! Denotes EOF condition
        INTEGER, INTENT(OUT)     :: IoStatus   ! I/O status code
    !
    ! !RETURN VALUE:
    !
        CHARACTER(LEN=MAXSTRLEN) :: Line       ! Single line from the input file

        !=================================================================
        ! Initialize
        !=================================================================
        EndOfFile = .FALSE.
        IoStatus  = 0
        Line      = ''

        !=================================================================
        ! Read data from the file
        !=================================================================

        ! Read a line from the file
        READ( fId, '(a)', IOSTAT=IoStatus ) Line

        ! IO Status < 0: EOF condition
        IF ( IoStatus < 0 ) THEN
           EndOfFile = .TRUE.
           RETURN
        ENDIF

        ! If desired, call StrSqueeze to strip leading and trailing blanks
        IF ( PRESENT( Squeeze ) ) THEN
           IF ( Squeeze ) THEN
              CALL StrSqueeze( Line )
           ENDIF
        ENDIF

      END FUNCTION ReadOneLine
    ! \brief Clean up a string
    !!
    !! \param  Str  Original string
    !! \return Cleaned-up string
    !!
    !! \ingroup core_modules
    !!!>
      FUNCTION CleanText( Str ) RESULT( CleanStr )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN) :: Str        ! Original string
    !
    ! !RETURN VALUE
    !
        CHARACTER(LEN=255)           :: CleanStr   ! Cleaned-up string

        ! Initialize
        CleanStr = Str

        ! Strip out non-printing characters (e.g. tabs)
        CALL CStrip    ( CleanStr           )

        ! Remove commas and quotes
        CALL StrRepl   ( CleanStr, ",", " " )
        CALL StrRepl   ( CleanStr, "'", " " )

        ! Remove leading and trailing spaces
        CALL StrSqueeze( CleanStr           )

      END FUNCTION CleanText
    ! \brief Print a string in multiple lines
    !!
    !! \param  Text        String to print
    !! \param  LineWidth   Width (characters) of lines
    !! \param  Delimiter   Delimiter between words
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE WordWrapPrint( Text, LineWidth, Delimiter )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*), INTENT(IN) :: Text        ! Text to print
        INTEGER,          INTENT(IN) :: LineWidth   ! Width (characters) of lines
        CHARACTER(LEN=1), OPTIONAL   :: Delimiter   ! Delimiter between words
    !
    ! !LOCAL VARIABLES:
    !
        ! Scalars
        INTEGER          :: C, S, B, Matches, Length

        ! Arrays
        INTEGER          :: BreakPts(100)
        INTEGER          :: SpaceLoc(500)

        ! Strings
        CHARACTER(LEN=1) :: Delim

        !=======================================================================
        ! WordWrapPrint begins here!
        !=======================================================================

        ! SpaceLoc is the array of where delimiters (usually the " "
        ! character) occur in the text, and S is its index.
        S           = 1
        SpaceLoc    = 0

        ! BreakPts is the array of where line breaks occur
        ! and B is its index.
        BreakPts    = 0
        B           = 1
        BreakPts(B) = 1

        ! Delimiter for separating words (will be the space character by default)
        IF ( PRESENT( Delimiter ) ) THEN
           Delim = Delimiter
        ELSE
           Delim = ' '
        ENDIF

        ! Find the Location of spaces in the text
        CALL CntMat( Text, ' ', Matches, SpaceLoc )

        ! Loop through the number of matches
        DO

           ! Move to the next delimiter location
           S = S + 1

           ! Compute the length of the line
           Length = SpaceLoc(S) - BreakPts(B)

           ! If the length of this segment is greater than the requested
           ! line length, store the position of this line break
           IF ( Length > LineWidth ) THEN
              B           = B             + 1
              BreakPts(B) = SpaceLoc(S-1) + 1
           ENDIF

           ! If we have exceeded the number of delimiters in the text, then set
           ! the last breakpoint at the end of the text and exit the loop.
           IF ( S > Matches ) THEN
              B           = B + 1
              BreakPts(B) = LEN_TRIM( Text ) + 1
              EXIT
           ENDIF

        ENDDO

        ! Print each line
        DO C = 1, B-1
           WRITE( 6, '(a)' ) Text( BreakPts(C):BreakPts(C+1)-1 )
        ENDDO

      END SUBROUTINE WordWrapPrint
    ! \brief Find unique elements in an array
    !!
    !! \param  vec  Input array
    !! \param vec_unique  Unique elements
    !!
    !! \ingroup core_modules
    !!!>
      SUBROUTINE Unique( vec, vec_unique )
    !
    ! !INPUT PARAMETERS:
    !
        CHARACTER(LEN=*),              INTENT(IN)    :: vec(:)
    !
    ! !OUTPUT PARAMETERS:
    !
        CHARACTER(LEN=*), ALLOCATABLE, INTENT(INOUT) :: vec_unique(:)

    !
    ! !LOCAL VARIABLES:
    !
        integer :: i,num
        logical,dimension(size(vec)) :: mask

        mask = .false.

        !=======================================================================
        ! Unique begins here!
        !=======================================================================

        ! Loop over all elements
        do i = 1, SIZE( vec )

           ! Assume that all valid array elements are located
           ! contiguously  Exit upon the encountering the
           ! first null character.
           IF ( LEN_TRIM( vec(I) ) == 0 ) EXIT

           !count the number of occurrences of this element:
           num = count( vec(i)==vec )

           if (num==1) then
              !there is only one, flag it:
              mask(i) = .true.
           else
              !flag this value only if it hasn't already been flagged:
              if (.not. any(vec(i)==vec .and. mask) ) mask(i) = .true.
           end if

        end do

        !return only flagged elements:
        IF ( ALLOCATED(  vec_unique ) ) DEALLOCATE( vec_unique )
        ALLOCATE( vec_unique(count(mask)) )
        vec_unique = PACK( vec, mask )

        !if you also need it sorted, then do so.
        ! For example, with slatec routine:
        !call ISORT (vec_unique, [0], size(vec_unique), 1)

      END SUBROUTINE Unique
    ! \brief Convert a character array to a string
    !!
    !! \param  charArray  Character array
    !! \param  N          Dimension of charArray
    !! \result     string     Output string
    !! \ingroup core_modules
    !!!>
      FUNCTION charArr2Str( charArray, N ) RESULT( string )
    !
    ! !INPUT PARAMETERS:
    !
        INTEGER,          INTENT(IN) :: N              ! Dimension of charArray
        CHARACTER(LEN=1), INTENT(IN) :: charArray(N)   ! Character array
    !
    ! !RETURN VALUE:
    !
        CHARACTER(LEN=N)             :: string         ! Output string
    !
    ! !LOCAL VARIABLES:
    !
        INTEGER :: C

        !=======================================================================
        ! charArray2str begins here!
        !=======================================================================

        ! Initialize the string
        string = ""


        ! Copy as much of the charArray to the string, until we hit the
        ! null byte (ASCII character 0), which denotes the end of characters
        DO C = 1, N
           IF ( chararray(C) == ACHAR(0) ) EXIT
           string(C:C) = chararray(C)
        ENDDO

      END FUNCTION charArr2Str
    ! \brief Convert a string to a character array
    !!
    !! \param  string  String
    !! \param  N       Length of string
    !! \result     charArray  Character array
    !! \ingroup core_modules
    !!!>
      FUNCTION str2CharArr( string, N ) RESULT( charArray )
    !
    ! !INPUT PARAMETERS:
    !
        INTEGER,          INTENT(IN) :: N              ! Length of string
        CHARACTER(LEN=N), INTENT(IN) :: string         ! Input string
    !
    ! !RETURN VALUE:
    !
        CHARACTER(LEN=1)             :: charArray(N)   ! Output character array
    !
    ! !LOCAL VARIABLES:
    !
       INTEGER :: C,  L

       ! Length of the string without trailing whitespace
       L = LEN_TRIM( string )

       ! Copy the non-whitespace characters to chararray
       DO C = 1, L
          chararray(C) = string(C:C)
       ENDDO

       ! Pad the remaining elements with the null byte
       charArray(L+1:) = ACHAR(0)

     END FUNCTION str2CharArr

    END MODULE CharPak_Mod
