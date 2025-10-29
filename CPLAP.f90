program cplap
!
!*********************************************************************
!
! This program reads in the chemical formula for a solid phase
! and its energy. It then asks for the limiting phases of the 
! combinations of the elements and their energies. These are 
! used to calculate the stability region in the space defined 
! by the chemical potentials of the constituent elements of the 
! solid. 
!
! The calculation is done by solving the system of linear equations
! resulting from the inequalities defined by the limiting constituent
! compound energies and the total energy expression of the solid. The
! chemical potential of the anionic species can be set as a free
! parameter.
!
! e.g. : 
!
! If the system is X Y2 Z with energy E ->
!
! we have mu_Z = E - mu_X - 2 mu_Y (*)
!
! with limiting phases X2 Z (energy E1)
!                      Y  Z (energy E2)
!
! we must have 2 mu_X + mu_Z < E1 (**)
!         and    mu_Y + mu_Z < E2 (**)
!
! in order that the system be grown (otherwise phases X2 Z or Y Z will 
! be formed)
!
! Also E < mu_X < 0 etc
!
! We can substitute in our expression (*) for mu_Z, then giving 3 
! equations with 2 unknowns. All combinations of the equations are 
! solved (if a solution exists for the combination) and a number of 
! intersection points are found (or not if the system has no solution).
! Those intersection points that satisfy the limiting inequalities (**)
! are the results we want, giving the corner points of the polyhedron
! defining the stability region in chemical potential space. (In the 
! example the polyhedron will be 2D).
!
! The limiting inequalities are returned, which can be used to 
! generate a plot in mathematica/matlab if the space is 3D or less.
! The corner points of the polyhedron defining the stability 
! region are also returned. The user can request a density of points
! (in units of energy), and a grid of points of that density in the 
! stability region are returned.
!
! The process can be started from a file. If this is done the user 
! is asked for a value of one of the chemical potentials (not the
! chemical potential that is a free parameter). The problem is then
! solved with this value of chemical potential, reducing the space 
! by one dimension.
!
! J. Buckeridge October 2011
!
! HACK - 10/03/2016: Check if solutions found are compatible with the
!        limiting inequalities has been moved to subroutine 
!        all_eqns_solve. Having it there avoids allocating very large 
!        arrays. A check has also been added a the program stops if too
!        much memory is used (i.e. if array 'results' has more than
!        nmax*nmax rows)
!
! BUGFIX - 29/10/2025: Problems identified with the restart (c) option,
!        where, if the reduced system is binary, there was an error in
!        some of the energies outputted. Fixed by subtracting the fixed
!        chemical potential value from the total energy. Also, minor 
!        problem with the neqns has been fixed, so that the correct lim-
!        iting inequalities are written. A test is required too to get
!        the correct energies for the inequalities related to the pure
!        elements.
!
!*********************************************************************
!
  implicit none
!
  integer, parameter :: nmax = 300
  real*8, parameter :: small = 1.d-12
  logical, parameter :: ltrue = .true.
  logical, parameter :: lfalse = .false.
!
  integer i, j, k, l, n, nspecies, nlimits, restart, num(nmax), nadd, nspecies_dep
  integer limittotno(nmax), limitnum(nmax,nmax), ndepend, neqns, numtmp
  integer numres, numsol, values(8), fixnum, str_num, comp_array(nmax), depnum
  integer input, extra, memtest
  integer nrel, func1, func2, rel_array(nmax), z_zero, cnt_z(nmax)
  real*8 mu(nmax), energy, limite(nmax), eqns(nmax,nmax)
  real*8 results(nmax*nmax,nmax), sum1, sum2, sum3, sum4
  real*8 tmp1, tmp2
  real*8 mumax(nmax), mumin(nmax), grid, gridpt(nmax), fixe
  real*8 bintemp, Ebin_min, Ebin_max, enertemp
  real*8 avgpos(nmax,3)
  character*1 read_tmp
  character*2 name(nmax), limitname(nmax,nmax), fixname, depname, tmp_str
  character*3 str2, str3
  character*4 str4
  character*5 zone, ineq_name(nmax), strdep
  character*8 date 
  character*9 ineq_str(nmax)
  character*10 time
  character*25 strtemp
  character*60 string, char
  character*100 str1, compounds(nmax)
  character*100 line
  character*200 eqn_str, formstr
  character*1000 lngstr
  logical check1, check2, check_eqns, check_low, check_ab, check_dep, check_lim
  logical check_fix, check_read
!
  open(unit=11, file="restart.dat", status="old", iostat=restart)
  open(unit=12, file="results.dat", status="replace")
!
! Write banner
!
  write(*,'(a)') "_____________________________________________________________"
  write(*,'(a)') "                                                             "
  write(*,'(a)') "     CCCC     PPPPP     LL         AAAA      PPPPP     !!    "
  write(*,'(a)') "    CC  CC    PP  PP    LL        AA  AA     PP  PP    !!    "
  write(*,'(a)') "    CC        PP  PP    LL        AA  AA     PP  PP    !!    "
  write(*,'(a)') "    CC        PPPPP     LL        AAAAAA     PPPPP     !!    "
  write(*,'(a)') "    CC  CC    PP        LL        AA  AA     PP              "
  write(*,'(a)') "     CCCC     PP        LLLLLL    AA  AA     PP        !!    "
  write(*,'(a)') "_____________________________________________________________"
  write(*,'(a)') "                                                             "
  write(*,'(a)') "         Chemical Potential Limits Analysis Program          "
  write(*,'(a)') "_____________________________________________________________"
  write(*,'(a)') "                                                             "
  write(*,'(a)') "  Version 1.1 2013                                           "
  write(*,'(a)') "  John Buckeridge (j.buckeridge@ucl.ac.uk)                   "
  write(*,'(a)') "_____________________________________________________________"
  check_ab = lfalse
  check_fix = lfalse
  if(restart /= 0) then
     write(*,*)
     write(*,'(a)') "************************   New run   ************************"
     write(*,'(a)') "_____________________________________________________________"
     write(*,*)
!
! Check for input file. If found read data from it
!
     open(unit=13, status="old", file="input.dat", iostat=input)
     if(input == 0) then
        write(*,'(a)') "Reading input from input.dat..."
        write(*,*)
!
! First read in the number of elements in the compound
!
        write(*,'(a)') "Reading in number of species in system..."
        do
           read(13,'(a100)') line
           if(line(1:1) /= "#") then
              read(line,*) nspecies
              exit
           endif
        enddo
        if(nspecies < 2) then
           write(*,'(a)') "SYNTAX ERROR reading in number of species (must be greater than one)"
           write(*,'(a)') "ERROR IS FATAL"
           goto 200
        endif
        write(*,'(a)') "...found ok"
        write(*,*)
!
! Now read in formula of compound
!
        write(*,'(a)') "Reading in formula of compound of interest..."
        do
           read(13,'(a100)') line
           if(line(1:1) /= "#") then
              read(line,*) (num(i),name(i),i=1,nspecies), energy
              exit
           endif
        enddo
!
! Check energy is negative
!
        if(energy > 0.d0) then
           write(*,'(a)') "FATAL ERROR: ENERGY OF SYSTEM IS POSITIVE!!"
           goto 200
        endif
!
! Check that element numbers in formula aren't negative
!
        do i=1,nspecies
           if(num(i) < 1) then
              write(*,'(a)') "FATAL ERROR: Incorrect syntax when entering numbers in chemical formula"
              write(*,'(a,i3)') "             Error occurred for element", i
              goto 200
           endif
        enddo
        write(*,'(a)') "...found ok"

     write(*,*)
     string = '(a7,   (1x,i2,1x,a2))'
     write(string(5:7),'(i3)') nspecies
     write(*,string) "System:", (num(i),name(i),i=1,nspecies)
     write(*,*)
!
! Read in which (if any) element is to be set as a dependent variable
!
        write(*,'(a)') "Reading in dependent variable..."
        do
           read(13,'(a100)') line
           if(line(1:1) /= "#") then
              read(line,*) depname
              exit
           endif
        enddo
!
! If none is set give a warning to the user
!
        if(depname(1:1) == "n") then
           write(*,'(a)') "WARNING: No depedent variable set. Solution may still be found -"
           write(*,'(a)') "         however only the intersection points with the compound"
           write(*,'(a)') "         are consistent with its formation. It is advisable to set"
           write(*,'(a)') "         a dependent variable."
           write(*,*)
           check_dep = lfalse
           nspecies_dep = nspecies
        else
!
! Check that the element is actually in the system
!
           check_dep = ltrue
           nspecies_dep = nspecies - 1
           j = 0
           do i=1,nspecies
              if(depname == name(i)) then
                 j = j+1
                 depnum = i
              endif
           enddo
           if(j == 0) then
              write(*,'(a)') "FATAL ERROR: Dependent variable element not found in system!"
              goto 200
           endif
!
! Change order of species in system so that dependent variable element is last
! in the list (ie element nspecies)
!
           tmp_str        = name(nspecies)
           name(nspecies) = depname
           name(depnum)   = tmp_str
           numtmp         = num(nspecies)
           num(nspecies)  = num(depnum)
           num(depnum)    = numtmp
           write(*,'(a)') "...found ok"
           write(*,*)
        endif
!
! Read in limiting phases and their energies. First read in the total number
!
        write(*,'(a)') "Reading in competing phases..."
        do
           read(13,'(a100)') line
           if(line(1:1) /= "#") then
              read(line,*) nlimits
              exit
           endif
        enddo
!
! Check that the number is not less than zero
!
        if(nlimits < 0) then
           write(*,'(a)') "FATAL ERROR: Found negative number of competing phases!"
           goto 200
        endif
!
! Check if there are no limiting compounds. If there are, read them in sequentially.
! Check that the number in the input file matches the stated total amount. If not warn
! user and reset total amount
!
        check_lim = ltrue
        if(nlimits == 0) then
           check_lim = lfalse
        else
           do i=1,nlimits
              check_read = lfalse
              do
                 read(13,'(a100)',end=5001) line
                 if(line(1:1) /= "#") then
                    read(line,*,end=5001) limittotno(i)
                    exit
                 endif
              enddo
              if(limittotno(i) > nspecies) then
                 write(*,'(a28,i3,a36)') "FATAL ERROR: Competing phase", i, "has more elements than the material!"
                 goto 200
              endif
              if(limittotno(i) < 1 ) then
                 write(*,'(a28,i3,a26)') "FATAL ERROR: Competing phase", i, "has less than one element!"
                 goto 200
              endif
              do
                 read(13,'(a100)',end=5001) line
                 if(line(1:1) /= "#") then
                    read(line,*,end=5001) (limitnum(i,j),limitname(i,j),j=1,limittotno(i)), limite(i)
                    exit
                 endif
              enddo
              check_read = ltrue
           enddo
5001       if(check_read .neqv. ltrue) then
              write(*,'(a)') "WARNING: Found less competing phases than supposed total amount"
              write(*,'(a29,i3)') "         Actual number found:", i-1
              nlimits = i-1
           endif
           write(*,'(a)') "...found ok"
           write(*,*)
!
! Check that elements in the compounds are compatible with those in the material
!
           write(*,'(a)') "Checking competing phase formulae..."
           do i=1,nlimits
              l = 0
              do j=1,limittotno(i)
                 do k=1,nspecies
                    if(limitname(i,j) == name(k)) l = l+1
                 enddo
              enddo
              if(l == 0) then
                 write(*,'(a)') "FATAL ERROR: Competing phase found with element not contained in &
                      &material!!"
                 write(*,'(a44,i3)') "             Error found in competing phase:", i
                 goto 200
              endif
           enddo
!
! Check that element numbers in compounds aren't negative
!
           do i=1,nlimits
              do j=1,limittotno(i)
                 if(limitnum(i,j) < 1) then
                    write(*,'(a)') "FATAL ERROR: Incorrect syntax when entering numbers in chemical formula"
                    write(*,'(a44,i3)') "             Error found in competing phase", i
                    goto 200
                 endif
              enddo
           enddo
           write(*,'(a)') "...all ok"
           write(*,*)
        endif
!
! If no input.dat file found then interactive mode is activated
!
     else
        write(*,'(a)') "No input.dat file found - beginning interactive mode..."
        write(*,*)
!
! Read in the number of species, chemical formula and energy
!    
        write(*,'(a)') "Enter number of elements in the material:"
70      read(*,*) nspecies
        if(nspecies < 2) then
           write(*,'(a)') "ERROR: System must consist of more than one element!!"
           write(*,'(a)') "Please re-enter number of elements in the material:"
           goto 70
        endif
!
        write(*,*)
        write(*,'(a)') "Enter chemical formula of material and its energy in the form"
        write(*,*)
        write(*,'(a)') "n_a A n_b B n_c C ... energy"
        write(*,*)
        write(*,'(a)') "where n_a etc. are the integer numbers for species A etc."
        write(*,'(a)') "(e.g. if the material is calcite then 1 Ca 1 C 3 O and then"
        write(*,'(a)') "the energy must be entered. Please keep the energy units"
        write(*,'(a)') "consistent when entering all data)"
        read(*,*) (num(i),name(i),i=1,nspecies), energy
        write(*,*)
!
! Check energy is negative
!
        if(energy > 0.d0) then
           write(*,'(a)') "FATAL ERROR: ENERGY OF SYSTEM IS POSITIVE!!"
           goto 200
        endif
!
! Check that element numbers in formula aren't negative
!
        do i=1,nspecies
           if(num(i) < 1) then
              write(*,'(a)') "FATAL ERROR: Incorrect syntax when entering numbers in chemical formula"
              write(*,'(a,i3)') "             Error occurred for element", i
              goto 200
           endif
        enddo
!
! Ask if an element is to be set as a dependent variable
!
        write(*,'(a)') "Do you wish to set an element as a dependent variable (y/n)?"
5000    read(*,'(a3)') str2
        if(str2(1:1) == "y" .or. str2(1:1) == "Y") then
           check_dep = ltrue
           nspecies_dep = nspecies - 1
           write(*,*)
           write(*,'(a)') "Enter element whose chemical potential will be the"
           write(*,'(a)') "dependent variable:"
8000       read(*,*) depname
           j = 0
           do i=1,nspecies
              if(depname == name(i)) then
                 j = j+1
                 depnum = i
              endif
           enddo
           if(j == 0) then
              write(*,'(a)') "ERROR: Element not found in system. Please enter another&
                   & element:"
              goto 8000
           endif
!
! Change order of species in system so that dependent variable element is last
! in the list (ie element nspecies)
!
           tmp_str        = name(nspecies)
           name(nspecies) = depname
           name(depnum)   = tmp_str
           numtmp         = num(nspecies)
           num(nspecies)  = num(depnum)
           num(depnum)    = numtmp
        elseif(str2(1:1) == "n" .or. str2(1:1) == "N") then
           write(*,*)
           write(*,'(a)') "WARNING: No depedent variable set. Solution may still be found -"
           write(*,'(a)') "         however only the intersection points with the compound"
           write(*,'(a)') "         are consistent with its formation. It is advisable to set"
           write(*,'(a)') "         a dependent variable."
           write(*,*)
           check_dep = lfalse
           nspecies_dep = nspecies
        else
           write(*,'(a)') "SYNTAX ERROR: please enter yes or no:"
           goto 5000
        endif
!
! Echo to user which chemical potential will be the dependent variable
!
        if(check_dep) then
           write(*,'(a5,1x,a)') "mu_"//trim(name(nspecies_dep + 1)), "will be the dependent &
                &variable"
           write(*,*)
        endif
!
! Read in limiting phases and their energies. Care will have to be taken when
! reading in each one, so the user will have to provide quite a lot of information
!
        write(*,'(a)') "---"
        write(*,'(a)')
        write(*,'(a)') "Now the competing phases and their energies must be entered."
        write(*,'(a)') "First enter the total number of competing phases:"
8001    read(*,*) nlimits
!
! Check that the number is not less than zero
!
        if(nlimits < 0) then
           write(*,'(a)') "ERROR: Number cannot be negative!"
           write(*,'(a)') "Enter number of competing phases again:"
           goto 8001
        endif
!
! Check if there are no limiting compounds. If so tell user that range of allowed
! chemical potentials is only limited by the allowed elemental values determined
! by the system energy
!
        check_lim = ltrue
        if(nlimits == 0) then
           write(*,'(a)') "No competing phases: element chemical potentials can take any"
           write(*,'(a)') "value in range between zero and system energy"
           check_lim = lfalse
        endif
        write(*,*)
        write(*,'(a)') "Now for each phase enter the number of elements, then the"
        write(*,'(a)') "formula and energy, in the same format as for the material"
        write(*,'(a)') "formula."
        write(*,*)
        do i=1,nlimits
           write(*,'(a48,i3)') "Enter number of elements in competing phase no. ", i
           read(*,*) limittotno(i)
!
! Check that the number of elements in the compound doesn't exceed the number
! of elements in the material, or is less than one
!
           if(limittotno(i) > nspecies) then
              write(*,'(a)') "ERROR: Competing phase has more elements than the material!!"
              write(*,*)
              write(*,'(a55,i3)') "Enter a new number of elements for competing phase no. ", i
              read(*,*) limittotno(i)
           endif
           if(limittotno(i) < 1 ) then
              write(*,'(a)') "ERROR: Competing phase cannot have less than one element!!"
              write(*,*)
              write(*,'(a55,i3)') "Enter a new number of elements for competing phase no. ", i
              read(*,*) limittotno(i)
           endif
           write(*,*)
           write(*,'(a57,i3)') "Enter chemical formula and energy of competing phase no. ", i
           read(*,*) (limitnum(i,j),limitname(i,j),j=1,limittotno(i)), limite(i)
           write(*,*)
        enddo
!
! Check that elements in the compounds are compatible with those in the material
!
        write(*,'(a)') "Checking competing phase formulae..."
        write(*,*)
        do i=1,nlimits
           l = 0
           do j=1,limittotno(i)
              do k=1,nspecies
                 if(limitname(i,j) == name(k)) l = l+1
              enddo
           enddo
           if(l == 0) then
              write(*,'(a)') "FATAL ERROR: Competing phase found with element not contained in &
                   &material!!"
              write(*,'(a44,i3)') "             Error found in competing phase ", i
              goto 200
           endif
        enddo
!
! Check that element numbers in compounds aren't negative
!
        do i=1,nlimits
           do j=1,limittotno(i)
              if(limitnum(i,j) < 1) then
                 write(*,'(a)') "FATAL ERROR: Incorrect syntax when entering numbers in chemical formula"
                 write(*,'(a44,i3)') "             Error found in competing phase ", i
                 goto 200
              endif
           enddo
        enddo
        write(*,'(a)') "...all ok"
        write(*,*)
     endif
!
! Construct matrix containing all linear equations. First go through all limiting
! compound formulae and determine the resulting linear equation for each
!
     eqns = 0.d0
     if(check_lim .eqv. lfalse) goto 100
     check1 = lfalse
     do i=1,nlimits
!
! First check if the dependent variable element is present. If so fill elements of 
! this line of the matrix appropriately, not yet taking into account the presence
! of other elements in the compound (skip this if no dependent variable set)
!
        if(check_dep) then
           do j=1,limittotno(i)
              if(limitname(i,j) == name(nspecies)) then
                 check1 = ltrue
                 ndepend = j
              endif
           enddo
        endif
        if(check1) then
           do j=1,nspecies_dep
              eqns(i,j) = eqns(i,j) - dble(limitnum(i,ndepend)) * (dble(num(j))/&
                   &dble(num(nspecies_dep + 1)))
           enddo
           eqns(i,nspecies_dep + 1) = limite(i) - dble(limitnum(i,ndepend)) * (energy/&
                &dble(num(nspecies_dep + 1)))
           check1 = lfalse
!
! If the dependent variable is not present the energy must still be put in the 
! last column
!
        else
           eqns(i,nspecies_dep + 1) = limite(i)
        endif
!
! Now fill matrix elements according to presence of other chemical elements besides 
! the dependent variable element. If the dependent variable element is present
! in the compound these numbers will be added to those determined previously. If not
! they will be added to zero (with matrix elements corresponding to chemical
! elements not present in the compound remaining as zero)
!
        do j=1,limittotno(i)
           do k=1,nspecies_dep
              if(limitname(i,j) == name(k)) then
                 eqns(i,k) = eqns(i,k) + dble(limitnum(i,j))
              endif
           enddo
        enddo
     enddo
!
! Now add equations corresponding to the limits imposed by the chemical potential
! ranges of the individual elements of the material (i.e. E < mu_X < 0). If no
! limiting compounds were inputted, definition of eqns array begins here
!
100  do i=1,nspecies_dep
        eqns(nlimits+i,i) = dble(num(i))
        eqns(nlimits+i+nspecies_dep,i) = dble(num(i))
        eqns(nlimits+i+nspecies_dep,nspecies_dep + 1) = energy
     enddo
!
! Now shift all equations forward one row in the matrix and add linear equation
! corresponding to the dependent variable equal to zero
!
     do i=nlimits+2*nspecies_dep,1,-1
        eqns(i+1,:) = eqns(i,:)
     enddo
     do i=1,nspecies_dep
        eqns(1,i) = dble(num(i))
     enddo
     eqns(1,nspecies_dep + 1) = energy
     neqns = nlimits + 2*nspecies_dep + 1
  else
!
! Restart file exists: read in data from restart file
!
     write(*,*)
     write(*,'(a)') "Found restart file..."
     write(*,*)
     do i=1,4
        read(11,*)
     enddo
     read(11,*) char, char, char, nspecies
     rewind(unit=11)
     read(11,*)
     read(11,*)
     read(11,*) char, char, char, (num(i),name(i),i=1,nspecies)
     write(*,*)
     string = '(a7,   (1x,i2,1x,a2))'
     write(string(5:7),'(i3)') nspecies
     write(*,string) "System:", (num(i),name(i),i=1,nspecies)
     read(11,*)
     read(11,*)
     read(11,*)
     read(11,*) char, char, char, char, nlimits
!
! Set check_lim
!
     if(nlimits > 0) then
        check_lim = ltrue
     else
        check_lim = lfalse
     endif
!
     do i=1,5
        read(11,*)
     enddo
     read(11,*) char, energy
     read(11,*)
     read(11,*) char, char, strdep
     if(strdep(1:1) == "n") then
        check_dep = lfalse
        nspecies_dep = nspecies
     else
        check_dep = ltrue
        nspecies_dep = nspecies - 1
        write(*,*)
        write(*,'(a19,2x,a5)') "Dependent variable:", strdep
     endif
     do i=1,5
        read(11,*)
     enddo
!
! Read in inequality relations to set up matrix eqns, containing the linear equations
! solved in the original run
!
     eqns = 0.d0
     do i=1,nlimits+1
        read(11,*) (eqns(i,j), char, j=1,nspecies_dep), char, eqns(i,nspecies_dep + 1)
     enddo
     do i=1,nspecies_dep
        read(11,*) eqns(nlimits+1+i,i), char, char, eqns(nlimits+1+i,nspecies_dep + 1)
     enddo
     do i=1,nspecies_dep
        read(11,*) eqns(nlimits+1+nspecies_dep+i,i), char, char, eqns(nlimits+1+&
             &nspecies_dep+i,nspecies_dep + 1)
     enddo
     read(11,*)
     read(11,*)
!
! Read in limiting compounds (if any)
!
     if(nlimits > 0) then
        limittotno = 0
        limitnum = 0
        limitname = " "
        limite = 0.d0
        do
           read(11,'(a1)') read_tmp
           if(read_tmp == "*") exit
        enddo
        read(11,*)
        read(11,*) char, char, numtmp, char, str1, char, char, char, limite(1)
!
! Parse out compound's elements from its name (contained in str1). Do this by 
! checking for lower case letters in the string, and for numbers
!
! Start with the first limiting compound
!
        str_num = len_trim(str1)
        i = 0
        do
           i = i+1
           check_low = lfalse
           if(i > str_num) exit
           if(i == str_num) then
              j = 1
           else
              j = 0
              do
                 j = j+1
                 if((iachar(str1(i+j:i+j)) >= iachar('A') .and. iachar(str1(i+j:i+j))&
                      <= iachar('Z')) .or. (j + i > str_num)) exit
              enddo
           endif
           if(j == 1) then
              limittotno(1) = limittotno(1) + 1
              limitname(1,limittotno(1)) = str1(i:i)
           else
              if(iachar(str1(i+1:i+1)) >= iachar('a') .and. iachar(str1(i+1:i+1)) &
                   <= iachar('z')) then
                 limittotno(1) = limittotno(1) + 1
                 limitname(1,limittotno(1)) = str1(i:i+1)
                 check_low = ltrue
              else
                 limittotno(1) = limittotno(1) + 1
                 limitname(1,limittotno(1)) = str1(i:i)
              endif
           endif
           if(check_low) then
              if(j > 2) then
                 do k=j,3,-1
                    limitnum(1,limittotno(1)) = limitnum(1,limittotno(1)) + (iachar(&
                         &str1(i+k-1:i+k-1)) - iachar('0')) * int(10.d0**(-k) * 10.d0**j)
                 enddo
              else
                 limitnum(1,limittotno(1)) = 1
              endif
           else
              if(j > 1) then
                 do k=j,2,-1
                    limitnum(1,limittotno(1)) = limitnum(1,limittotno(1)) + (iachar(&
                         &str1(i+k-1:i+k-1)) - iachar('0')) * int(10.d0**(-k) * 10.d0**j)
                 enddo
              else
                 limitnum(1,limittotno(1)) = 1
              endif
           endif
           i = i + j - 1
        enddo
!
! Now the other limiting compounds (if any)
!
        if(nlimits > 1) then
           do
              do
                 read(11,'(a1)') read_tmp
                 if(read_tmp == "_") exit
              enddo
              str1 = " "
              read(11,*)
              read(11,*) char, char, numtmp, char, str1, char, char, char, limite&
                   &(numtmp)
              str_num = len_trim(str1)
              i = 0
              do
                 i = i+1
                 check_low = lfalse
                 if(i > str_num) exit
                 if(i == str_num) then
                    j = 1
                 else
                    j = 0
                    do
                       j = j+1
                       if((iachar(str1(i+j:i+j)) >= iachar('A') .and. iachar(str1(i+j&
                            &:i+j)) <= iachar('Z')) .or. (j + i > str_num)) exit
                    enddo
                 endif
                 if(j == 1) then
                    limittotno(numtmp) = limittotno(numtmp) + 1
                    limitname(numtmp,limittotno(numtmp)) = str1(i:i)
                 else
                    if(iachar(str1(i+1:i+1)) >= iachar('a') .and. iachar(str1(i+1:i+1))&
                         <= iachar('z')) then
                       limittotno(numtmp) = limittotno(numtmp) + 1
                       limitname(numtmp,limittotno(numtmp)) = str1(i:i+1)
                       check_low = ltrue
                    else
                       limittotno(numtmp) = limittotno(numtmp) + 1
                       limitname(numtmp,limittotno(numtmp)) = str1(i:i)
                    endif
                 endif
                 if(check_low) then
                    if(j > 2) then
                       do k=j,3,-1
                          limitnum(numtmp,limittotno(numtmp)) = limitnum(numtmp,&
                               &limittotno(numtmp)) + (iachar(str1(i+k-1:i+k-1)) -&
                               &iachar('0')) * int(10.d0**(-k) * 10.d0**j)
                       enddo
                    else
                       limitnum(numtmp,limittotno(numtmp)) = 1
                    endif
                 else
                    if(j > 1) then
                       do k=j,2,-1
                          limitnum(numtmp,limittotno(numtmp)) = limitnum(numtmp,&
                               &limittotno(numtmp)) + (iachar(str1(i+k-1:i+k-1)) -&
                               &iachar('0')) * int(10.d0**(-k) * 10.d0**j)
                       enddo
                    else
                       limitnum(numtmp,limittotno(numtmp)) = 1
                    endif
                 endif
                 i = i + j - 1
              enddo
              if(numtmp == nlimits) exit
           enddo
        endif
     endif
!
! Three restart options are available. Ask the user which is appropriate
!
     write(*,*)
     write(*,*)
     write(*,'(a)') "From restart file can:"
     write(*,*)
     write(*,'(a)') "     (a) Change which element's chemical potential is the"
     write(*,'(a)') "         dependent variable (or set a dependent variable if"
     write(*,'(a)') "         not done previously)"
     write(*,'(a)') "     (b) Enter additional competing phases (can be read from"
     write(*,'(a)') "         a file extra_phases.dat)"
     write(*,'(a)') "     (c) Set the value of chemical potential for a particular"
     write(*,'(a)') "         element to a constant"
     write(*,*)
     write(*,'(a)') "Enter appropriate choice (a/b/c):"
90   read(*,'(a1)') read_tmp
!
     if(read_tmp(1:1) == 'a' .or. read_tmp(1:1) == 'A') then
!
! Dependent variable must be changed. Ask user which element's chemical potential
! will be the new dependent variable. Also output file will be slightly different
! to the (b) option case, so need to set a check here
!
        check_ab = ltrue
        check_dep = ltrue
        nspecies_dep = nspecies - 1
        write(*,*)
        write(*,'(a)') "Enter element whose chemical potential will be the (new)"
        write(*,'(a)') "dependent variable:"
80      read(*,*) fixname
        j = 0
        do i=1,nspecies
           if(fixname == name(i)) then
              j = j+1
              fixnum = i
           endif
        enddo
        if(j == 0) then
           write(*,'(a)') "ERROR: Element not found in system. Please enter another&
                & element:"
           goto 80
        endif
!
! Change order of species in system so that new dependent variable element is last
! in the list (ie element nspecies)
!
        tmp_str        = name(nspecies)
        name(nspecies) = fixname
        name(fixnum)   = tmp_str
        numtmp         = num(nspecies)
        num(nspecies)  = num(fixnum)
        num(fixnum)    = numtmp
!
! Now create eqns array from the limiting compounds and system, as was done for the
! original calculation
!
        eqns = 0.d0
        if(check_lim .eqv. lfalse) goto 110
        check1 = lfalse
        do i=1,nlimits
!
! First check if the dependent variable element is present. If so fill elements of 
! this line of the matrix appropriately, not yet taking into account the presence
! of other elements in the compound
!
           do j=1,limittotno(i)
              if(limitname(i,j) == name(nspecies_dep + 1)) then
                 check1 = ltrue
                 ndepend = j
              endif
           enddo
           if(check1) then
              do j=1,nspecies_dep
                 eqns(i,j) = eqns(i,j) - dble(limitnum(i,ndepend)) * (dble(num(j))/&
                      &dble(num(nspecies_dep + 1)))
              enddo
              eqns(i,nspecies_dep + 1) = limite(i) - dble(limitnum(i,ndepend)) * (energy/&
                   &dble(num(nspecies_dep + 1)))
              check1 = lfalse
!
! If the dependent variable is not present the energy must still be put in the 
! last column
!
           else
              eqns(i,nspecies_dep + 1) = limite(i)
           endif
!
! Now fill matrix elements according to presence of other chemical elements besides 
! the dependent variable element. If the dependent variable element is present
! in the compound these numbers will be added to those determined previously. If not
! they will be added to zero (with matrix elements corresponding to chemical
! elements not present in the compound remaining as zero)
!
           do j=1,limittotno(i)
              do k=1,nspecies_dep
                 if(limitname(i,j) == name(k)) then
                    eqns(i,k) = eqns(i,k) + dble(limitnum(i,j))
                 endif
              enddo
           enddo
        enddo
!
! Now add equations corresponding to the limits imposed by the chemical potential
! ranges of the individual elements of the material (i.e. E < mu_X < 0). If no
! limiting compounds were inputted, definition of eqns array begins here
!
110     do i=1,nspecies_dep
           eqns(nlimits+i,i) = dble(num(i))
           eqns(nlimits+i+nspecies_dep,i) = dble(num(i))
           eqns(nlimits+i+nspecies_dep,nspecies_dep + 1) = energy
        enddo
!
! Now shift all equations forward one row in the matrix and add linear equation
! corresponding to the dependent variable equal to zero
!
        do i=nlimits+2*nspecies_dep,1,-1
           eqns(i+1,:) = eqns(i,:)
        enddo
        do i=1,nspecies_dep
           eqns(1,i) = dble(num(i))
        enddo
        eqns(1,nspecies_dep + 1) = energy
        neqns = nlimits + 2*nspecies_dep + 1
!
! If option b is selected, new limiting compounds must be read in
!
     elseif(read_tmp(1:1) == 'b' .or. read_tmp(1:1) == 'B') then
!
        check_ab = ltrue
!
! Check for file containing additional phases. If it exists read data from it
!
        open(unit=14, status="old", file="extra_phases.dat", iostat=extra)
        if(extra == 0) then
           write(*,'(a)') "Reading additional competing phases from extra_phases.dat..."
           write(*,*)
           write(*,'(a)') "Reading total number of additional phases..."
           do
              read(14,'(a100)') line
              if(line(1:1) /= "#") then
                 read(line,*) nadd
                 exit
              endif
           enddo
!
! Check that the number is not less than zero
!
           if(nadd < 0) then
              write(*,'(a)') "FATAL ERROR: Found negative number of additional competing phases!"
              goto 200
           endif
           write(*,'(a)') "...found ok"
           write(*,*)
!
! Check if there are no limiting compounds. If there are, read them in sequentially.
! Check that the number in the input file matches the stated total amount. If not warn
! user and reset total amount
!
           if(nadd == 0) then
              write(*,'(a)') "Zero additional competing phases found"
              goto 130
           else
              write(*,'(a)') "Reading in formulae of each addtional competing phase..."
              write(*,*)
              do i=nlimits+1,nlimits+nadd
                 check_read = lfalse
                 do
                    read(14,'(a100)',end=5002) line
                    if(line(1:1) /= "#") then
                       read(line,*,end=5002) limittotno(i)
                       exit
                    endif
                 enddo
                 if(limittotno(i) > nspecies) then
                    write(*,'(a28,i3,a36)') "FATAL ERROR: Competing phase", i, "has more elements than the material!"
                    goto 200
                 endif
                 if(limittotno(i) < 1 ) then
                    write(*,'(a28,i3,a26)') "FATAL ERROR: Competing phase", i, "has less than one element!"
                    goto 200
                 endif
                 do
                    read(14,'(a100)',end=5002) line
                    if(line(1:1) /= "#") then
                       read(line,*,end=5002) (limitnum(i,j),limitname(i,j),j=1,limittotno(i)), limite(i)
                       exit
                    endif
                 enddo
                 check_read = ltrue
              enddo
5002          if(check_read .neqv. ltrue) then
                 write(*,'(a)') "WARNING: Found less additional competing phases than supposed total amount"
                 write(*,'(a29,i3)') "         Actual number found:", i-1
                 nadd = i-1
              endif
              write(*,'(a)') "...found ok"
              write(*,*)
!
! Check that elements in the compounds are compatible with those in the material
!
              write(*,'(a)') "Checking additional competing phase formulae..."
              do i=nlimits+1,nlimits+nadd
                 l = 0
                 do j=1,limittotno(i)
                    do k=1,nspecies
                       if(limitname(i,j) == name(k)) l = l+1
                    enddo
                 enddo
                 if(l == 0) then
                    write(*,'(a)') "FATAL ERROR: Competing phase found with element not contained in &
                         &material!!"
                    write(*,'(a44,i3)') "             Error found in competing phase:", i
                    goto 200
                 endif
              enddo
!
! Check that element numbers in compounds aren't negative
!
              do i=nlimits+1,nlimits+nadd
                 do j=1,limittotno(i)
                    if(limitnum(i,j) < 1) then
                       write(*,'(a)') "FATAL ERROR: Incorrect syntax when entering numbers in chemical formula"
                       write(*,'(a44,i3)') "             Error found in competing phase", i
                       goto 200
                    endif
                 enddo
              enddo
              write(*,'(a)') "...all ok"
              write(*,*)
           endif
        else
!
! File not found. Enter values interactively
!
           write(*,'(a)') "extra_phases.dat file not found - entering interactive mode..."
           write(*,*)
           write(*,'(a)') "Enter number of additional competing phases:"
120        read(*,*) nadd
           if(nadd < 0 .or. nadd > 1000) then
              write(*,*)
              write(*,'(a)') "ERROR: number cannot be negative!!"
              write(*,'(a)') "Please re-enter number:"
              goto 120
           endif
!
           if(nadd == 0) then
              write(*,'(a)') "No additional competing phases"
              goto 130
           endif
           write(*,*)
           write(*,'(a)') "For each competing phase enter the number of elements, then the"
           write(*,'(a)') "formula and energy, in the same format as for the material"
           write(*,'(a)') "formula."
           write(*,*)
           do i=1,nadd
              write(*,'(a59,i3)') "Enter number of elements in additional competing phase no. ", i
              read(*,*) limittotno(nlimits+i)
!
! Check that the number of elements in the compound doesn't exceed the number
! of elements in the material, or is less than 1
!
              if(limittotno(nlimits+i) > nspecies) then
                 write(*,'(a)') "ERROR: Competing phase has more elements than the material!!"
                 write(*,*)
                 write(*,'(a66,i3)') "Enter a new number of elements for additional &
                      &competing phase no. ", i
                 read(*,*) limittotno(nlimits+i)
              endif
              if(limittotno(nlimits+i) < 1 ) then
                 write(*,'(a)') "ERROR: Competing phase cannot have less than one element!!"
                 write(*,*)
                 write(*,'(a55,i3)') "Enter a new number of elements for competing phase no. ", i
                 read(*,*) limittotno(nlimits+i)
              endif
              write(*,*)
              write(*,'(a67,i3)') "Enter chemical formula and energy of additional &
                   &competing phase no. ", i
              read(*,*) (limitnum(nlimits+i,j),limitname(nlimits+i,j),j=1,limittotno&
                   &(nlimits+i)), limite(nlimits+i)
              write(*,*)
           enddo
!
! Check that elements in the compounds are compatible with those in the material
!
           write(*,'(a)') "Checking competing phase formulae..."
           write(*,*)
           do i=1,nadd
              l = 0
              do j=1,limittotno(nlimits+i)
                 do k=1,nspecies
                    if(limitname(nlimits+i,j) == name(k)) l = l+1
                 enddo
              enddo
              if(l == 0) then
                 write(*,'(a)') "FATAL ERROR: Competing phase found with element not contained in &
                      &material!!"
                 write(*,'(a44,i3)') "             Error found in competing phase ", i
                 goto 200
              endif
           enddo
!
! Check that element numbers in compounds aren't negative
!
           do i=1,nadd
              do j=1,limittotno(nlimits+i)
                 if(limitnum(nlimits+i,j) < 1) then
                    write(*,'(a)') "FATAL ERROR: Incorrect syntax when entering numbers in chemical formula"
                    write(*,'(a44,i3)') "             Error found in competing phase ", i
                    goto 200
                 endif
              enddo
           enddo
           write(*,'(a)') "...all ok"
           write(*,*)
!
! Close if testing if extra_phases.dat exists
!
        endif
!
! Reset nlimits value and check_lim
!
        nlimits = nlimits + nadd
        check_lim = ltrue
!
! Now create eqns array from the limiting compounds and system, as was done for the
! original calculation
!
130     eqns = 0.d0
        if(check_lim .eqv. lfalse) goto 140
        check1 = lfalse
        do i=1,nlimits
!
! First check if the dependent variable element is present. If so fill elements of 
! this line of the matrix appropriately, not yet taking into account the presence
! of other elements in the compound (skip this if no dependent variable set)
!
           if(check_dep) then
              do j=1,limittotno(i)
                 if(limitname(i,j) == name(nspecies)) then
                    check1 = ltrue
                    ndepend = j
                 endif
              enddo
           endif
           if(check1) then
              do j=1,nspecies_dep
                 eqns(i,j) = eqns(i,j) - dble(limitnum(i,ndepend)) * (dble(num(j))/&
                      &dble(num(nspecies_dep + 1)))
              enddo
              eqns(i,nspecies_dep + 1) = limite(i) - dble(limitnum(i,ndepend)) * (energy/&
                   &dble(num(nspecies_dep + 1)))
              check1 = lfalse
!
! If the dependent variable is not present the energy must still be put in the 
! last column
!
           else
              eqns(i,nspecies_dep + 1) = limite(i)
           endif
!
! Now fill matrix elements according to presence of other chemical elements besides 
! the dependent variable element. If the dependent variable element is present
! in the compound these numbers will be added to those determined previously. If not
! they will be added to zero (with matrix elements corresponding to chemical
! elements not present in the compound remaining as zero)
!
           do j=1,limittotno(i)
              do k=1,nspecies_dep
                 if(limitname(i,j) == name(k)) then
                    eqns(i,k) = eqns(i,k) + dble(limitnum(i,j))
                 endif
              enddo
           enddo
        enddo
!
! Now add equations corresponding to the limits imposed by the chemical potential
! ranges of the individual elements of the material (i.e. E < mu_X < 0). If no
! limiting compounds were inputted, definition of eqns array begins here
!
140     do i=1,nspecies_dep
           eqns(nlimits+i,i) = dble(num(i))
           eqns(nlimits+i+nspecies_dep,i) = dble(num(i))
           eqns(nlimits+i+nspecies_dep,nspecies_dep + 1) = energy
        enddo
!
! Now shift all equations forward one row in the matrix and add linear equation
! corresponding to the dependent variable equal to zero
!
        do i=nlimits+2*nspecies_dep,1,-1
           eqns(i+1,:) = eqns(i,:)
        enddo
        do i=1,nspecies_dep
           eqns(1,i) = dble(num(i))
        enddo
        eqns(1,nspecies_dep + 1) = energy
        neqns = nlimits + 2*nspecies_dep + 1
!
! Otherwise option c is requested
!
     elseif(read_tmp(1:1) == 'c' .or. read_tmp(1:1) == 'C') then   
!
! If the system is binary this option is not available
!
        if(nspecies <= 2) then
           write(*,*)
           write(*,'(a)') "ERROR: binary system found - option (c) not available!!"
           write(*,'(a)') "       (solution is trivial)"
           write(*,*)
           write(*,'(a)') "Please select option a or b:"
           goto 90
        endif
!
! Request from user element for which chemical potential will be fixed
!
        write(*,*)
        write(*,'(a)') "Enter element for which chemical potential will be fixed:"
60      read(*,*) fixname
        j = 0
        do i=1,nspecies
           if(fixname == name(i)) then
              j = j+1
              if(check_dep .and. i == nspecies) then
                 write(*,'(a)') "ERROR: The dependent variable cannot be held fixed.&
                      &Please"
                 write(*,'(a)') "enter another element:"
                 goto 60
              else
                 fixnum = i
              endif
           endif
        enddo
        if(j == 0) then
           write(*,'(a)') "ERROR: Element not found in system. Please enter another &
                &element:"
           goto 60
        endif
        check_fix = ltrue
!
! Request value at which this element's chemical potential should be set
!
        write(*,*)
        write(*,'(a)') "Enter fixed value of chemical potential for this element:"
        read(*,*) fixe
!
! Re-stack eqns matrix as now one of the chemical potentials has a fixed value
!
! BUGFIX 29/10/2025: "+ 1" was missing from the following line. Needed to get the last 
! limiting inequality printed
!
        neqns = nlimits + 2*nspecies_dep + 1
        do i=1,neqns
! 
! BUGFIX 29/10/2025: Added if statement for cases where we are looking at inequalities
! involving the pure elements. In these cases, the "eqns(i,fixnum)" is zero, but the energy
! still needs to be changed by the fixed chemical potential
!
           if (i > nlimits + nspecies + 1) then 
              eqns(i,nspecies_dep+1) = eqns(i,nspecies_dep+1) - num(fixnum) * fixe
           else         
              eqns(i,nspecies_dep+1) = eqns(i,nspecies_dep+1) - eqns(i,fixnum) * fixe
           endif
        enddo
        do i=fixnum,nspecies_dep
           eqns(:,i) = eqns(:,i+1)
        enddo
        do i=nlimits+1+fixnum,neqns
           eqns(i,:) = eqns(i+1,:)
        enddo
        do i=nlimits+nspecies_dep+fixnum,neqns
           eqns(i,:) = eqns(i+1,:)
        enddo
!
! Arrays 'name' and 'num' must also be changed
!
! BUGFIX 29/10/2025: energy needs to have the fixed chemical potential value subtracted from it
! (of course, multiplied by the appropriate factor, the number related to the fixed chemical potential 
! element)
!
        energy = energy - num(fixnum) * fixe
!
        do i=fixnum,nspecies_dep
           num(i) = num(i+1)
           name(i) = name(i+1)
        enddo
        nspecies = nspecies - 1
        nspecies_dep = nspecies_dep - 1
        neqns = neqns - 2
!
! Else incorrect input from user. Ask for correct input
!
     else
        write(*,'(a)') "SYNTAX ERROR: please enter a, b or c:"
        goto 90
     endif
!
! Close if-else checking for restart file
!
  endif
!
! Write to file the limiting relations, including a header with system and date and
! time, energy and dependent variable
!
  call date_and_time(date,time,zone,values)
  call date_and_time(DATE=date,ZONE=zone)
  call date_and_time(TIME=time)
  call date_and_time(VALUES=values)
  write(12,'(a)') "**************************************************************&
       &********************"
  write(12,*)
  if(restart /= 0 .or. check_ab) then
     string = '(a19,   (1x,i2,1x,a2))'
     write(string(6:8),'(i3)') nspecies
     write(12,string) "Results for system:", (num(i),name(i),i=1,nspecies)
  else
     string = '(a24,   (1x,i2,1x,a2),3x,a1,1x,a5,2x,a1,2x,f9.4,1x,a1)'
     write(string(6:8),'(i3)') nspecies
     write(12,string) "Results- reduced system:", (num(i),name(i),i=1,nspecies)&
          &,"(", "mu_"//fixname, "=", fixe, ")"
  endif
  write(12,*)
  write(12,'(a18,i3)') "Number of species:", nspecies
  write(12,*)
  write(12,'(a27,i3)') "Number of competing phases:", nlimits
  write(12,*)
  write(12,'(a5,2x,I2,a1,I2,a1,I4,3x,I2,a1,I2,a1,I2)') 'DATE:', &
       &values(3), '/', values(2), '/', values(1), values(5), ':', values(6), ':',&
       &values(7)
  write(12,*) 
  write(12,'(a)') "**************************************************************&
       &********************"
  write(12,*) 
  write(12,'(a7,e21.13,1x,a2)') "Energy:", energy, "eV"
  write(12,*)
  if(check_dep) then
     write(12,'(a19,2x,a5)') "Dependent variable:", "mu_"//name(nspecies)
  else
     write(12,'(a19,2x,a5)') "Dependent variable:", "none"
  endif
  write(12,*)
  write(12,'(a)') "**************************************************************&
       &********************"
  write(12,*)
  write(12,'(a)') "Limiting inequalities:"
  write(12,*)
  do i=1,nlimits+1
     do j=1,nspecies_dep
        ineq_name(j) = "mu_"//name(j)
        write(ineq_str(j)(1:9),'(f9.4)') eqns(i,j)
     enddo
     if(i == 1) then
        string = '(   (a8,1x,a5,2x),a1,2x,f9.4)'
        write(string(2:4),'(i3)') nspecies_dep
        if(check_dep) then
           write(12,string) (ineq_str(j),ineq_name(j),j=1,nspecies_dep), ">", &
                &  eqns(i,nspecies_dep + 1)
        else
           write(12,string) (ineq_str(j),ineq_name(j),j=1,nspecies_dep), "=", &
                &  eqns(i,nspecies_dep + 1)
        endif
     else
        string = '(   (a8,1x,a5,2x),a1,2x,f9.4)'
        write(string(2:4),'(i3)') nspecies_dep
        write(12,string) (ineq_str(j),ineq_name(j),j=1,nspecies_dep), "<", &
             &eqns(i,nspecies_dep + 1)
     endif
  enddo
  do i=nlimits+2,neqns
     do j=1,nspecies_dep
        if(eqns(i,j) /= 0.d0) then
           ineq_name(j) = "mu_"//name(j)
           write(ineq_str(j)(1:9),'(f9.4)') eqns(i,j)
        else
           ineq_name(j) = "     "
           ineq_str(j) = "        "
        endif
     enddo
     if(i > nlimits + nspecies_dep + 1) then
        string = '(   (a8,1x,a5,2x),a1,2x,f9.4)'
        write(string(2:4),'(i3)') nspecies_dep
        write(12,string) (ineq_str(j),ineq_name(j),j=1,nspecies_dep), ">", &
           &  eqns(i,nspecies_dep + 1)
     else
        string = '(   (a8,1x,a5,2x),a1,2x,f9.4)'
        write(string(2:4),'(i3)') nspecies_dep
        write(12,string) (ineq_str(j),ineq_name(j),j=1,nspecies_dep), "<", &
             &eqns(i,nspecies_dep + 1)
     endif
  enddo
     write(12,*)
     write(12,'(a)') "**************************************************************&
          &********************"
!
! If system is binary there is only one independent variable, and results will
! have to be determined in a different way to the more general case (there will
! be no systems of equations to solve). If no limiting compounds were supplied 
! the solution is trivial
!
  if(nspecies_dep == 1) then
!
     write(*,*)
     write(*,'(a)') "Determining solution for binary system..."
     write(*,*)
!
     if(nlimits == 0) then
        Ebin_min = energy / dble(num(1))
        Ebin_max = 0.d0
     else
!
! There is a solution if the limits on the independent variable are consistent.
! Redefine the limits on the independent variable according to the limiting 
! compounds, and check if there is a solution
!
        Ebin_min = energy / dble(num(1))
        Ebin_max = 0.d0
        do i=2,nlimits+1
           bintemp = eqns(i,nspecies) / eqns(i,nspecies-1)
           if(eqns(i,nspecies-1) < 0.d0) then
              if(bintemp > Ebin_min) Ebin_min = bintemp
           else
              if(bintemp < Ebin_max) Ebin_max = bintemp
           endif
        enddo
     endif
     if(Ebin_min > Ebin_max) then
!
! No solution exists
!
        write(*,*)
        write(*,'(a)') "******************************************************"
        write(*,'(a)') "**                                                  **"
        write(*,'(a)') "**        NO SOLUTIONS FOUND FOR MATERIAL!!!        **"
        write(*,'(a)') "**                                                  **"
        write(*,'(a)') "**             MATERIAL IS UNSTABLE!!!              **"
        write(*,'(a)') "**                                                  **"
        write(*,'(a)') "******************************************************"
        write(*,*)
        write(12,*)
        write(12,'(a)') "System is unstable - no solutions exist"
        write(12,*)
        write(12,'(a)') "**************************************************************&
             &********************"   
     else
        numsol = 2
        results(1,1) = Ebin_min
        results(2,1) = Ebin_max
!
! Write solution to file
!
        write(12,*)
        write(*,'(a)') "Solution found - material is thermodynamically stable!"
        write(*,*)
        write(12,'(a)') "Solution for binary system:"
        write(12,*)
        write(12,'(2x,f9.4,2x,a1,2x,a5,2x,a1,2x,f9.4)') Ebin_min, "<", "mu_"//name(1), "<"&
             &, Ebin_max
        write(12,*)
        write(12,*)
        write(12,'(2x,a5,2x,a1,2x,f9.4,4x,a2,4x,a5,2x,a2,2x,f9.4)') "mu_"//name(1), "=", &
             &Ebin_min, "->", "mu_"//name(2), "=", (energy - dble(num(1)) * Ebin_min) / &
             &dble(num(2))
        write(12,*)
        write(12,'(2x,a5,2x,a1,2x,f9.4,4x,a2,4x,a5,2x,a2,2x,f9.4)') "mu_"//name(1), "=", &
             &Ebin_max, "->", "mu_"//name(2), "=", (energy - dble(num(1)) * Ebin_max) / &
             &dble(num(2))
        write(12,*)
        write(12,*)
        write(12,'(a)') "Equation for system:"
        write(12,*)
        write(12,'(2x,f9.4,2x,a5,2x,a1,2x,f9.4,2x,a5,2x,a1,2x,f9.4)') dble(num(1)), "mu_"&
             &//name(1), "+", dble(num(2)), "mu_"//name(2), "=", energy
        write(12,*)
        write(12,'(a)') "**************************************************************&
             &********************"
        write(*,'(a)') "...done"
     endif
     if(nlimits > 0) then
        goto 160
     else
        goto 170
     endif
  endif
!
! Call routine to solve all combinations of linear equations for systems with more 
! than two elements
!
  numsol = 0
  call all_eqns_solve(eqns(1:neqns,1:nspecies_dep + 1),neqns,nspecies_dep + 1,&
       &results(:,1:nspecies_dep),nlimits,numsol)
!
! If no solutions are found tell the user
!
  if(numsol == 0) then
     write(*,*)
     write(*,'(a)') "******************************************************"
     write(*,'(a)') "**                                                  **"
     write(*,'(a)') "**        NO SOLUTIONS FOUND FOR MATERIAL!!!        **"
     write(*,'(a)') "**                                                  **"
     write(*,'(a)') "**             MATERIAL IS UNSTABLE!!!              **"
     write(*,'(a)') "**                                                  **"
     write(*,'(a)') "******************************************************"
     write(*,*)
     write(12,*)
     write(12,'(a)') "System is unstable - no solutions exist"
     write(12,*)
     write(12,'(a)') "**************************************************************&
          &********************"   
     goto 200
!
! Otherwise write the solutions to file
!
  else
     write(*,'(a)') "Solution found - material is thermodynamically stable!"
     write(*,*)
!
! First check for redundant solutions and remove from list before writing to file
!
     do i=1,numsol
        j = i
        do
           j = j+1
           if(j > numsol) exit
           l = 0
           do k=1,nspecies_dep
              if(results(i,k) == results(j,k)) l = l+1
           enddo
           if(l == nspecies_dep) then
              do n=j,numsol-1
                 results(n,:) = results(n+1,:)
              enddo
              numsol = numsol - 1
              j = j-1
           endif
        enddo
     enddo
!
! Now write solutions
!
     write(12,*)
     write(12,'(a)') "Intersection points in chemical potential space:"
     write(12,*)
     string = '(   (3x,a5,3x),a1,3x,a5,1x,a1)'
     write(string(2:4),'(i3)') nspecies-1
     write(12,string) ("mu_"//name(i),i=1,nspecies-1), "|", "mu_"//name(nspecies), "="
     str1 = "-----------"
     str1 = trim(str1)
     do i=1,nspecies-2
        str1 = "-----------"//str1
     enddo
     write(str1((nspecies-1)*11+1:(nspecies-1)*11+14),'(a13)') "|------------"
     string = '(a   )'
     write(str2(1:3),'(i3)') (nspecies * 11) + 1 
     str3 = adjustl(str2)
     write(string(3:5),'(a3)') str3
     write(12,string) str1
     do i=1,numsol
        string = '(   (1x,f9.4,1x),a1,1x,f9.4)'
        write(string(2:4),'(i3)') nspecies-1
        sum1 = 0.d0
        do j=1,nspecies-1
           sum1 = sum1 + dble(num(j)) * results(i,j)
        enddo
        enertemp = (energy - sum1) / dble(num(nspecies))
        write(12,string) (results(i,j),j=1,nspecies-1), "|", enertemp 
     enddo
     write(12,*)
     write(12,'(a)') "**************************************************************&
          &********************"   
  endif
  if(nlimits > 0) then
     goto 160
  else
     goto 170
  endif
!
! Now write equations to file corresponding to each limiting compound, including
! the corresponding intersection points if any (do this even if there is no solution
! for the system) 
!
  if(nlimits == 0) goto 150
160 do i=1,nlimits
     write(12,*)
!
! Create string with correct chemical name for limiting compound
!
     str_num = 0
     str1 = " "
     do j=1,limittotno(i)
        if(len_trim(limitname(i,j)) > 1) then
           write(str1(str_num+1:str_num+2),'(a2)') limitname(i,j)
        else
           write(str1(str_num+1:str_num+1),'(a1)') limitname(i,j)(1:1)
        endif
        str_num = len_trim(str1)
        if(limitnum(i,j) > 1) then
           if(limitnum(i,j) < 10) then
              write(str1(str_num+1:str_num+1),'(i1)') limitnum(i,j)
           elseif(limitnum(i,j) > 9 .and. limitnum(i,j) < 100) then
              write(str1(str_num+1:str_num+2),'(i2)') limitnum(i,j)
           else
              write(str1(str_num+1:str_num+3),'(i3)') limitnum(i,j)
           endif
        endif
        str_num = len_trim(str1)
     enddo
     compounds(i+1) = str1
!
! Write compound name as a header, with its energy
!
     string = "(a15,1x,i3,1x,a1,3x,a   ,3x,a10,1x,f9.4,1x,a1)"
     write(str2,'(i3)') str_num
     str2 = adjustl(str2)
     write(string(22:24),'(a3)') str2
     write(12,string) "Competing phase", i, ":", str1, "( energy =", limite(i), ")"
     write(12,'(a23)') "-----------------------"
     write(12,*)
!
! Create string consisting of relevant equation for compound
!
     eqn_str = " "
     check_eqns = ltrue
     if(eqns(i+1,1) /= 0.d0) then
        write(eqn_str(1:17),'(f9.4,1x,a5,2x)') eqns(i+1,1), "mu_"//name(1)
        check_eqns = lfalse
     endif
     do j=2,nspecies_dep
        if(check_eqns) then
           if(eqns(i+1,j) /= 0.d0) then
              write(eqn_str(1:17),'(f9.4,1x,a5,2x)') eqns(i+1,j), "mu_"//name(j)
              check_eqns = lfalse
           endif
        else
           if(eqns(i+1,j) /= 0.d0) then
              str_num = len_trim(eqn_str)
              if(eqns(i+1,j) < 0.d0) then
                 write(eqn_str(str_num+3:str_num+23),'(a1,2x,f9.4,1x,a5,2x)') "-", &
                      &-eqns(i+1,j), "mu_"//name(j)
              else
                 write(eqn_str(str_num+3:str_num+23),'(a1,2x,f9.4,1x,a5,2x)') "+", &
                      &eqns(i+1,j), "mu_"//name(j)
              endif
              check_eqns = lfalse
           endif
        endif
     enddo
     str_num = len_trim(eqn_str)
     write(eqn_str(str_num+3:str_num+14),'(a1,2x,f9.4)') "=", eqns(i+1,nspecies_dep + 1)
!
! Write the equation to file
!
     write(12,'(a)') "Equation:"
     write(12,*)
     str_num = len_trim(eqn_str)
     string = "(2x,a   )"
     write(str2,'(i3)') str_num
     str2 = adjustl(str2)
     write(string(6:8),'(a3)') str2
     write(12,string) eqn_str
     write(12,*)
!
! Determine which intersection points are relevant for each compound and write them
! to file
!
     write(12,'(a)') "Relevant intersection points:"
     write(12,*)
     comp_array = 0
     l = 0
     if(numsol > 0) then
        do j=1,numsol
           check2 = lfalse
           sum1 = 0.d0
           do k=1,nspecies_dep
              sum1 = sum1 + results(j,k) * eqns(i+1,k)
           enddo
           if(abs(sum1 - eqns(i+1,nspecies_dep + 1)) < small) then
              check2 = ltrue
           endif
!
! Count how many there are and put their order in an array
!
           if(check2) then
              l = l+1
              comp_array(l) = j
           endif
        enddo
     endif
!
! Now write to file. If there are none then say so
!
     if(l == 0 .or. numsol == 0) then
        write(12,'(a)') "None"
     else
        string = '(   (3x,a5,3x),a1,3x,a5,1x,a2)'
        write(string(2:4),'(i3)') nspecies-1
        write(12,string) ("mu_"//name(j),j=1,nspecies-1), "|", "mu_"//name(nspecies)&
             &, ">="
        str1 = "-----------"
        str1 = trim(str1)
        do j=1,nspecies-2
           str1 = "-----------"//str1
        enddo
        write(str1((nspecies-1)*11+1:(nspecies-1)*11+14),'(a13)') "|------------"
        string = '(a   )'
        write(str2(1:3),'(i3)') (nspecies * 11) + 1 
        str3 = adjustl(str2)
        write(string(3:5),'(a3)') str3
        write(12,string) str1
        do j=1,l
           string = '(   (1x,f9.4,1x),a1,1x,f9.4)'
           write(string(2:4),'(i3)') nspecies-1
           sum1 = 0.d0
           do k=1,nspecies-1
              sum1 = sum1 + dble(num(k)) * results(comp_array(j),k)
           enddo
           enertemp = (energy - sum1) / dble(num(nspecies))
           write(12,string) (results(comp_array(j),k),k=1,nspecies-1), "|", enertemp 
        enddo
     endif
     write(12,'(a)') "______________________________________________________________&
          &____________________"
  enddo

!
! Repeat for the limits imposed by the system itself
!
170 write(12,*)
!
! Create string with correct chemical name for limiting compound
!
  str1 = " "
  str_num = 0
  do j=1,nspecies
     if(len_trim(name(j)) > 1) then
        write(str1(str_num+1:str_num+2),'(a2)') name(j)
     else
        write(str1(str_num+1:str_num+1),'(a1)') name(j)(1:1)
     endif
     str_num = len_trim(str1)
     if(num(j) > 1) then
        if(num(j) < 10) then
           write(str1(str_num+1:str_num+1),'(i1)') num(j)
        elseif(num(j) > 9 .and. num(j) < 100) then
           write(str1(str_num+1:str_num+2),'(i2)') num(j)
        else
           write(str1(str_num+1:str_num+3),'(i3)') num(j)
        endif
     endif
     str_num = len_trim(str1)
  enddo
  compounds(1) = str1
!
! Write compound name as a header, with its energy
!
  string = "(a23,3x,a   ,3x,a10,1x,f9.4,1x,a1)"
  write(str2,'(i3)') str_num
  str2 = adjustl(str2)
  write(string(10:12),'(a3)') str2
  write(12,string) "Limits from system    :", str1, "( energy =", energy, ")"
  write(12,'(a23)') "-----------------------"
  write(12,*)
!
! Create string consisting of relevant equation for compound
!
  eqn_str = " "
  check_eqns = ltrue
  if(eqns(1,1) /= 0.d0) then
     write(eqn_str(1:17),'(f9.4,1x,a5,2x)') eqns(1,1), "mu_"//name(1)
     check_eqns = lfalse
  endif
  do j=2,nspecies_dep
     if(check_eqns) then
        if(eqns(1,j) /= 0.d0) then
           write(eqn_str(1:17),'(f9.4,1x,a5,2x)') eqns(1,1), "mu_"//name(j)
           check_eqns = lfalse
        endif
     else
        if(eqns(1,j) /= 0.d0) then
           str_num = len_trim(eqn_str)
           if(eqns(1,j) < 0.d0) then
              write(eqn_str(str_num+3:str_num+23),'(a1,2x,f9.4,1x,a5,2x)') "-", &
                   &-eqns(1,j), "mu_"//name(j)
           else
              write(eqn_str(str_num+3:str_num+23),'(a1,2x,f9.4,1x,a5,2x)') "+", &
                   &eqns(1,j), "mu_"//name(j)
           endif
           check_eqns = lfalse
        endif
     endif
  enddo
  str_num = len_trim(eqn_str)
  write(eqn_str(str_num+3:str_num+14),'(a1,2x,f9.4)') "=", eqns(1,nspecies_dep + 1)
!
! Write the equation to file
!
  write(12,'(a)') "Equation:"
  write(12,*)
  str_num = len_trim(eqn_str)
  string = "(2x,a   )"
  write(str2,'(i3)') str_num
  str2 = adjustl(str2)
  write(string(6:8),'(a3)') str2
  write(12,string) eqn_str
  write(12,*)
!
! Determine which intersection points are relevant for the system and write them
! to file
!
  write(12,'(a)') "Relevant intersection points:"
  write(12,*)
  comp_array = 0
  l = 0
  if(numsol > 0) then
     do j=1,numsol
        check2 = lfalse
        sum1 = 0.d0
        do k=1,nspecies_dep
           sum1 = sum1 + results(j,k) * eqns(1,k)
        enddo
        if(abs(sum1 - eqns(1,nspecies_dep + 1)) < small) then
           check2 = ltrue
        endif
!
! Count how many there are and put their order in an array
!
        if(check2) then
           l = l+1
           comp_array(l) = j
        endif
     enddo
  endif
!
! Now write to file. If there are none then say so
!
  if(l == 0 .or. numsol == 0) then
     write(12,'(a)') "None"
  else
     string = '(   (3x,a5,3x),a1,3x,a5,1x,a2)'
     write(string(2:4),'(i3)') nspecies-1
     write(12,string) ("mu_"//name(j),j=1,nspecies-1), "|", "mu_"//name(nspecies)&
          &, ">="
     str1 = "-----------"
     str1 = trim(str1)
     do j=1,nspecies-2
        str1 = "-----------"//str1
     enddo
     write(str1((nspecies-1)*11+1:(nspecies-1)*11+14),'(a13)') "|------------"
     string = '(a   )'
     write(str2(1:3),'(i3)') (nspecies * 11) + 1 
     str3 = adjustl(str2)
     write(string(3:5),'(a3)') str3
     write(12,string) str1
     do j=1,l
        string = '(   (1x,f9.4,1x),a1,1x,f9.4)'
        write(string(2:4),'(i3)') nspecies-1
        sum1 = 0.d0
        do k=1,nspecies-1
           sum1 = sum1 + dble(num(k)) * results(comp_array(j),k)
        enddo
        enertemp = (energy - sum1) / dble(num(nspecies))
        write(12,string) (results(comp_array(j),k),k=1,nspecies-1), "|", enertemp 
     enddo
  endif
  write(12,*)
  write(12,'(a)') "**************************************************************&
       &********************"
150 continue
!
! Ask user if grid of points in polyhedron is required. Skip this if the system is
! binary
!
  if(nspecies <= 2 .or. numsol == 0) goto 180
  write(*,*)
  write(*,'(a)') "Print grid of points enclosed within stability polyhedron (y/n)?"
50 read(*,'(a3)') str2
  if(str2(1:1) == "y" .or. str2(1:1) == "Y") then
     open(unit=15, file="grid.dat", status="replace")
     write(*,*)
     write(*,'(a)') "Enter grid dimension in units of energy (e.g. 0.1, &
          &0.01 etc.):"
     read(*,*) grid
!
! Check that grid dimension is not too small
!
10   continue
     if(grid < small*1.d5) then
        write(*,'(a)') "ERROR: Grid dimension is too small. Please enter a larger &
             &value:"
        read(*,*) grid
        goto 10
     endif
!
! Create grid. First create regular grid (rectangular or cubic or higher dimensional
! depending on nspecies-1) using the max and min values of each mu_X. Then check each
! point in this grid to see if it is compatible with the limiting inequalities. If it
! is write it to file
!
     write(*,*)
     write(*,'(a)') "Determining grid and writing to file..."
     write(*,*)
     mumax = -1.d8
     mumin =  1.d8
     do i=1,nspecies-1
        do j=1,numsol
           if(results(j,i) > mumax(i)) then
              mumax(i) = results(j,i)
           endif
           if(results(j,i) < mumin(i)) then
              mumin(i) = results(j,i)
           endif
        enddo
     enddo
!
! Write to file the header and the polyhedron corner points
!
     call date_and_time(date,time,zone,values)
     call date_and_time(DATE=date,ZONE=zone)
     call date_and_time(TIME=time)
     call date_and_time(VALUES=values)
     write(15,'(a)') "**************************************************************&
          &********************"
     write(15,*)
     if(restart /= 0) then
        string = '(a19,   (1x,i2,1x,a2))'
        write(string(6:8),'(i3)') nspecies
        write(15,string) "Results for system:", (num(i),name(i),i=1,nspecies)
     else
        string = '(a24,   (1x,i2,1x,a2),3x,a1,1x,a5,2x,a1,2x,f9.4,1x,a1)'
        write(string(6:8),'(i3)') nspecies
        write(15,string) "Results- reduced system:", (num(i),name(i),i=1,nspecies),&
             &"(", "mu_"//fixname, "=", fixe, ")"
     endif
     write(15,*)
     write(15,'(a5,2x,I2,a1,I2,a1,I4,3x,I2,a1,I2,a1,I2)') 'DATE:', &
          values(3), '/', values(2), '/', values(1), values(5), ':', values(6), ':',&
          &values(7)
     write(15,*) 
     write(15,'(a9,i3,1x,a39)') "The first", numsol, "points are the polyhedron corner&
         & points"
     write(15,*)
     write(15,'(a)') "**************************************************************&
          &********************" 
     write(15,*)
     string = '(   (3x,a5,3x),a1,3x,a5,1x,a2)'
     write(string(2:4),'(i3)') nspecies-1
     write(15,string) ("mu_"//name(i),i=1,nspecies-1), "|", "mu_"//name(nspecies), ">="
     str1 = "-----------"
     str1 = trim(str1)
     do i=1,nspecies-2
        str1 = "-----------"//str1
     enddo
     write(str1((nspecies-1)*11+1:(nspecies-1)*11+14),'(a13)') "|------------"
     string = '(a   )'
     write(str2(1:3),'(i3)') (nspecies * 11) + 1 
     str3 = adjustl(str2)
     write(string(3:5),'(a3)') str3
     write(15,string) str1
     do i=1,numsol
        string = '(   (1x,f9.4,1x),a1,1x,f9.4)'
        write(string(2:4),'(i3)') nspecies-1
        sum1 = 0.d0
        do j=1,nspecies-1
           sum1 = sum1 + dble(num(j)) * results(i,j)
        enddo
        enertemp = (energy - sum1) / dble(num(nspecies))
        write(15,string) (results(i,j),j=1,nspecies-1), "|", enertemp 
     enddo
!
! Start at the min value for each mu_X, increment by grid, and test if the resulting
! point is within the polyhedron
!
     do i=1,nspecies-1
        gridpt(i) = mumin(i)
     enddo
20   continue
     check2 = ltrue
     do i=1,neqns
        sum1 = 0.d0
        do j=1,nspecies-1
           sum1 = sum1 + gridpt(j) * eqns(i,j)
        enddo
        if(i == 1 .or. i > nlimits + nspecies) then
           if(sum1 - eqns(i,nspecies) < -small) then
              check2 = lfalse
           endif
        else
           if(sum1 - eqns(i,nspecies) > small) then
              check2 = lfalse
           endif
        endif
     enddo
!
! Write compatible grid points to file
!
     if(check2) then
        string = '(   (1x,f9.4,1x),a1,1x,f9.4)'
        write(string(2:4),'(i3)') nspecies-1
        sum1 = 0.d0
        do j=1,nspecies-1
           sum1 = sum1 + dble(num(j)) * gridpt(j)
        enddo
        enertemp = (energy - sum1) / dble(num(nspecies))
        write(15,string) (gridpt(j),j=1,nspecies-1), "|", enertemp
     endif
     i = nspecies-1
30   gridpt(i) = gridpt(i) + grid
     if(gridpt(i) < mumax(i)) then
        goto 20
     else
        gridpt(i) = mumin(i)
        i = i-1
        if(i < 1) then
           goto 40
        else
           goto 30
        endif
     endif
40   write(*,'(a)') "...done"
     close(unit=15)
!
  elseif(str2(1:1) == "n" .or. str2(1:1) == "N") then
     write(*,*)
     write(*,'(a)') "No grid printed"
     write(*,*)
  else
     write(*,'(a)') "SYNTAX ERROR: please enter yes or no:"
     goto 50
  endif
!
! If System of independent variables is 2D or 3D a plot file can be generated
! for mathematica. For the 3D case a plot file is generated only if there is a 
! solution
! A plotting file for gnuplot is also produced in both cases
!
180 if(nspecies_dep == 3) then
     open(unit=16,status='replace',file="3Dplot.dat")
     write(*,*)
     write(*,'(a)') "Chemical potential space is 3D. Writing plotting files..." 
     write(16,'(a20,f7.3,a2,2(a1,f7.3,a2),a1,f7.3,a2)') "plot1:=RegionPlot3D[", &
          &eqns(1,1), name(1), ("+", eqns(1,i), name(i), i=2,3), ">", eqns(1,nspecies_dep + 1),&
          "&&"
     do i=2,nlimits
        write(16,'(20x,f7.3,a2,2(a1,f7.3,a2),a1,f7.3,a2)') eqns(i,1), name(1), ("+", &
             &eqns(i,j), name(j), j=2,3), "<", eqns(i,nspecies_dep + 1),"&&"
     enddo
     write(16,'(20x,f7.3,a2,2(a1,f7.3,a2),a1,f7.3,a1)') eqns(nlimits+1,1), name(1), &
          &("+", eqns(nlimits+1,j), name(j), j=2,3), "<", eqns(nlimits+1,nspecies_dep + 1),","
     write(16,'(20x,3(a1,a2,a1,f7.3,a4))') ("{", name(i), ",", energy/dble(num(i)), ",0},"&
          &,i=1,3)
     write(16,'(20x,a62,a11)') "PlotRange->Full,PlotPoints->100,PlotStyle->Directive[Yellow],"&
          &, "Mesh->None,"
     write(16,'(20x,a12,2(a2,a1),a2,a47)') "AxesLabel->{", (name(i),",",i=1,2), name(3), &
          "},Ticks->Automatic,ViewPoint->{Pi,3*Pi/2,Pi/3}]"
     write(16,'(a20,f7.3,a9,f7.3,a9,f7.3,a14)') "plot2:=ListPlot3D[{{", energy/dble(num(1))&
          &, ",0,0},{0,", energy/dble(num(2)), ",0},{0,0,", energy/dble(num(3)), &
          &"}},Mesh->None,"
     write(16,'(18x,a35)') "PlotStyle->Directive[Opacity[0.1]]]"
     write(16,'(a17)') "Show[plot1,plot2]"
     close(unit=16)
!
! Now write gnuplot file
!
     open(unit=16,status='replace',file="3Dplot.plt")
     write(16,'(a)') "set view equal xyz"
     write(16,'(a)') "set view 60, 220"
     write(16,'(a)') "set border 4095"
     write(16,'(a)') "set xyplane 0"
     write(16,*)
!
! Go through each relevant compound, including the main compound. Check if any solution points
! correspond to that compound. If so, it must be included in the plot. Find the intersection 
! lines with all other relevant compounds, and convert to conditional expressions on x and y.
! Surfaces are defined with reference to z. The midpoints of the surfaces are used to work out
! the conditional expressions. Surfaces parallel to z cannot be plotted - instead their 
! midpoint is plotted
!
! Define numerical values for the functions so that they can be plotted. Also need to keep track
! of number of relevant compounds and their order (in an array)
!
     func1 = 1000
     func2 = 100
     nrel = 0
     rel_array = 0
!
! Go through each compound
!
     do i=1,nlimits+1
!
! Check if any solutions correspond to that compound, and put their order in an array if they
! do. Also find the mid-point of the solution points and store that
!
        comp_array = 0
        l = 0
        sum2 = 0.d0
        sum3 = 0.d0
        sum4 = 0.d0
        do j=1,numsol
           check2 = lfalse
           sum1 = 0.d0
           do k=1,nspecies_dep
              sum1 = sum1 + results(j,k) * eqns(i,k)
           enddo
           if(abs(sum1 - eqns(i,nspecies_dep + 1)) < small) then
              check2 = ltrue
           endif
           if(check2) then
              l = l+1
              comp_array(l) = j
              sum2 = sum2 + results(j,1)
              sum3 = sum3 + results(j,2)
              sum4 = sum4 + results(j,3)
           endif
        enddo
!
! If there are solutions for this compound, count it and put its order in an array, and 
! work out the midpoint
!
        if(l > 0) then
           nrel = nrel + 1
           rel_array(nrel) = i
           avgpos(nrel,1) = sum2/real(l)
           avgpos(nrel,2) = sum3/real(l)
           avgpos(nrel,3) = sum4/real(l)
        endif
     enddo
!
! Go through each compound that has relevant intersection points. For each case, write 
! the formula of the corresponding plane in terms of z, and go through all other relevant
! compounds and find the intersection lines. Put these in a conditional function on x and
! y. The less than or greater than sign is determined from the position of the mid-point.
! Planes parallel to z (i.e. with zero z coefficient) must be treated differently - 
! for these, count how many there are and find their order. They will be added to the splot
! command differently so that their midpoint only will be plotted 
!
     z_zero = 0
     cnt_z = 0
     do i=1,nrel
        func1 = func1 + 1
        func2 = func2 + 1
        lngstr = ""
!
! Skip planes parallel to z, but count them and store their order
!
        if(eqns(rel_array(i),3) == 0.d0) then
           z_zero = z_zero+1
           cnt_z(z_zero) = i
           cycle
        endif
!
! First write the equation of the surface for the compound
!
        write(16,'(a1,i4,a7,f9.4,a3,f9.4,a6,f9.4,a4)') "f", func1, "(x,y)=(", &
             &eqns(rel_array(i),nspecies_dep+1)/eqns(rel_array(i),3), " + ", &
             &-eqns(rel_array(i),1)/eqns(rel_array(i),3), "* x + ", &
             &-eqns(rel_array(i),2)/eqns(rel_array(i),3), "* y)" 
!
! Now create string to write the conditional expressions
!
        write(lngstr,'(a1,i3,a7)') "f", func2, "(x,y)=("
!
! Write the first condition separately to the rest as it has a slightly different format
!
        k = 0
        do
           k = k+1
           if(k /= i) then
              if(eqns(rel_array(k),3) == 0.d0) then
                 sum1 = 0.d0
                 sum2 = 0.d0
                 do j=1,2
                    sum1 = sum1 + eqns(rel_array(k),j) * avgpos(i,j)
                 enddo
                 sum2 = eqns(rel_array(k),nspecies_dep+1)
                 if(sum1 <= sum2) then
                    str4 = "*y<="
                 else
                    str4 = "*y>="
                 endif
                 write(string,'(f9.4,a5,f9.4,a4,f9.4)') eqns(rel_array(k),1),&
                      &"*x + ", eqns(rel_array(k),2), str4, sum2
                 lngstr = trim(lngstr)//trim(string)
                 exit
              else
                 sum1 = 0.d0
                 sum2 = 0.d0
                 do j=1,2
                    sum1 = sum1 + ((eqns(rel_array(k),j)/eqns(rel_array(k)&
                         &,3))-(eqns(rel_array(i),j)/eqns(rel_array(i),3))) * avgpos(i,j)
                 enddo
                 sum2 = (eqns(rel_array(k),nspecies_dep+1)/eqns(rel_array(k),3))-(eqns(&
                      &rel_array(i),nspecies_dep+1)/eqns(rel_array(i),3))
                 if(sum1 <= sum2) then
                    str4 = "*y<="
                 else
                    str4 = "*y>="
                 endif
                 write(string,'(f9.4,a5,f9.4,a4,f9.4)') (eqns(rel_array(k),1)/eqns(rel_array(k)&
                      &,3))-(eqns(rel_array(i),1)/eqns(rel_array(i),3)), "*x + ", (eqns(rel_array&
                      &(k),2)/eqns(rel_array(k),3))-(eqns(rel_array(i),2)/eqns(rel_array(i),3)),&
                      &str4, sum2
                 lngstr = trim(lngstr)//trim(string)
                 exit
              endif
           endif
!
! If this conditional loop goes wrong, give an escape and report an error
!
           if(k > nrel) then
              write(*,'(a)') "ERROR: problem writing 3D gnuplot plotting file"
              close(unit=16)
              goto 200
           endif
        enddo
!
! Now write the rest of the conditional expressions to the string, with "&&" in front of each one
!
        if(k < nrel) then
           do j=k+1,nrel
              if(j == i) cycle
              if(eqns(rel_array(j),3) == 0.d0) then
                 sum1 = 0.d0
                 sum2 = 0.d0
                 do l=1,2
                    sum1 = sum1 + eqns(rel_array(j),l) * avgpos(i,l)
                 enddo
                 sum2 = eqns(rel_array(j),nspecies_dep+1)
                 if(sum1 <= sum2) then
                    str4 = "*y<="
                 else
                    str4 = "*y>="
                 endif
                 write(string,'(a2,f9.4,a5,f9.4,a4,f9.4)') "&&", eqns(rel_array(j),1),&
                      &"*x + ", eqns(rel_array(j),2), str4, sum2
                 lngstr = trim(lngstr)//trim(string)
              else
                 sum1 = 0.d0
                 sum2 = 0.d0
                 do l=1,2
                    sum1 = sum1 + ((eqns(rel_array(j),l)/eqns(rel_array(j)&
                         &,3))-(eqns(rel_array(i),l)/eqns(rel_array(i),3))) * avgpos(i,l)
                 enddo
                 sum2 = (eqns(rel_array(j),nspecies_dep+1)/eqns(rel_array(j),3))-(eqns(&
                      &rel_array(i),nspecies_dep+1)/eqns(rel_array(i),3))
                 if(sum1 <= sum2) then
                    str4 = "*y<="
                 else
                    str4 = "*y>="
                 endif
                 write(string,'(a2,f9.4,a5,f9.4,a4,f9.4)') "&&", (eqns(rel_array(j),1)/eqns(rel_array&
                      &(j),3))-(eqns(rel_array(i),1)/eqns(rel_array(i),3)), "*x + ", (eqns(rel_array&
                      &(j),2)/eqns(rel_array(j),3))-(eqns(rel_array(i),2)/eqns(rel_array(i),3)), &
                      &str4, sum2
                 lngstr = trim(lngstr)//trim(string)
              endif
           enddo
        endif
!
! Complete the string and write it to file
!
        write(string,'(a3,i4,a9)') ")?f", func1, "(x,y):1/0"
        lngstr = trim(lngstr)//trim(string)
        write(16,'(a)') trim(lngstr)
     enddo
!
! More settings for the plot
!
     write(16,*)
     write(16,'(a)') "set isosample 300, 300"
     write(16,'(a)') "set hidden3d offset 0"
     write(16,'(a)') "set key outside left center"
     write(16,'(a12,f9.4,a3)') "set xrange [", eqns(1,nspecies_dep + 1), ":0]"
     write(16,'(a12,f9.4,a3)') "set yrange [", eqns(1,nspecies_dep + 1), ":0]"
     write(16,'(a12,f9.4,a3)') "set zrange [", eqns(1,nspecies_dep + 1), ":0]"
     write(16,'(a12,a2,a1)') "set xlabel '", name(1), "'"
     write(16,'(a12,a2,a1)') "set ylabel '", name(2), "'"
     write(16,'(a12,a2,a1)') "set zlabel '", name(3), "'"
     write(16,*)
!
! Set arrows to mark the limits imposed by the system
!
     write(16,'(a26,f9.4,a6,f9.4,a12)') "set arrow nohead from 0,0,", eqns(1,nspecies_dep+1)&
          &/eqns(1,3), " to 0,", eqns(1,nspecies_dep+1)/eqns(1,2), ",0 lt 3 lw 2"
     write(16,'(a26,f9.4,a4,f9.4,a14)') "set arrow nohead from 0,0,", eqns(1,nspecies_dep+1)&
          &/eqns(1,3), " to ", eqns(1,nspecies_dep+1)/eqns(1,1), ",0,0 lt 3 lw 2"
     write(16,'(a24,f9.4,a6,f9.4,a14)') "set arrow nohead from 0,", eqns(1,nspecies_dep+1)/&
          &eqns(1,2), ",0 to ", eqns(1,nspecies_dep+1)/eqns(1,1), ",0,0 lt 3 lw 2"
     write(16,*)
!
! Now write splot command, making sure formatting is correct. Write it as a string first, as this 
! allows each case to be looped over more easily. For cases with zero z coefficient,
! a point will be plotted instead of a surface
!

     lngstr = ""
     write(lngstr,'(a19,a,a2)') "splot f101(x,y) t '", trim(compounds(rel_array(1))), "' "
     do i=2,nrel
        if(z_zero /= 0) then
           check2 = lfalse
           do j=1,z_zero
              if(i == cnt_z(j)) check2 = ltrue
           enddo
           if(check2) cycle
        endif
        write(string,'(a3,i3,a9,a,a1)') ", f", 100+i, "(x,y) t '", trim(compounds(rel_array(i))), "'"
        lngstr = trim(lngstr)//trim(string)
     enddo
     if(z_zero /= 0) then
        do i=1,z_zero
           write(string,'(a13,a,a1)') ", '-' w p t '", trim(compounds(rel_array(cnt_z(i)))), "'"
           lngstr = trim(lngstr)//trim(string)
        enddo
     endif
!
! Write splot string to file
!
     write(16,'(a)') trim(lngstr)
!
! If there any planes parallel to z, write their midpoints to be plotted
!
     if(z_zero /= 0) then
        do i=1,z_zero
           write(16,'(3(f9.4,1x))') (avgpos(cnt_z(i),j),j=1,3)
           write(16,'(a1)') "e"
        enddo
     endif
!
! Finish
!
     close(unit=16)
     write(*,*)
     write(*,'(a)') "...done -> 3Dplot.plt can be loaded in gnuplot"
     write(*,'(a)') "        -> 3Dplot.dat can be pasted into Mathematica"
!
! If system is 2D prepare plotting files
!
  elseif(nspecies_dep == 2) then
     open(unit=17,status='replace',file="2Dplot.dat")
     write(*,*)
     write(*,'(a)') "Chemical potential space is 2D. Writing files for plotting..."
     write(17,'(a18,f7.3,a2,a1,f7.3,a2,a1,f7.3,a2)') "plot1:=RegionPlot[", &
          &eqns(1,1), name(1), "+", eqns(1,2), name(2), ">", eqns(1,nspecies_dep + 1),&
          "&&"
     do i=2,nlimits
        write(17,'(18x,f7.3,a2,a1,f7.3,a2,a1,f7.3,a2)') eqns(i,1), name(1), "+", &
             &eqns(i,2), name(2), "<", eqns(i,nspecies_dep + 1),"&&"
     enddo
     write(17,'(18x,f7.3,a2,a1,f7.3,a2,a1,f7.3,a1)') eqns(nlimits+1,1), name(1), &
          &"+", eqns(nlimits+1,2), name(2), "<", eqns(nlimits+1,nspecies_dep + 1),","
     write(17,'(18x,2(a1,a2,a1,f7.3,a4))') ("{", name(i), ",", energy/dble(num(i)), ",0},"&
          &,i=1,2)
     write(17,'(18x,a32,a11)') "PlotRange->Full,PlotPoints->100,"
     write(17,'(18x,a12,a2,a1,a2,a2)') "AxesLabel->{", name(1),",", name(2), "}]"
     write(17,'(a14,f7.3,a1,f7.3,a2,a2)') "plot2:=Plot[{{", eqns(1,nspecies_dep + 1)/eqns(1,2), "-", &
          &eqns(1,1)/eqns(1,2), name(1), "},"
     do i=1,nlimits-1
        write(17,'(14x,a2,f7.3,a1,f7.3,a2,a2)') "{", eqns(i+1,nspecies_dep + 1)/eqns(i+1,2), &
             &"-", eqns(i+1,1)/eqns(i+1,2), name(1), "},"
     enddo
     write(17,'(14x,a2,f7.3,a1,f7.3,a2,a3)') "{", eqns(nlimits+1,nspecies_dep + 1)/eqns(nlimits+1&
          &,2), "-", eqns(nlimits+1,1)/eqns(nlimits+1,2), name(1), "}},"
     write(17,'(14x,a1,a2,a1,f7.3,a4)') "{", name(1), ",", eqns(1,nspecies_dep + 1)/eqns(1,1), &
          &",0}]"
     write(17,'(a17)') "Show[plot1,plot2]"
     close(unit=17)
!
! Now write gnuplot file
!
     open(unit=17,status='replace',file="2Dplot.plt")
     write(17,'(a)') "set multiplot"
     write(17,*)
     write(17,'(a)') "set lmargin at screen 0.225"
     write(17,'(a)') "set rmargin at screen 0.775"
     write(17,'(a)') "set tmargin at screen 0.875"
     write(17,'(a)') "set bmargin at screen 0.125"
     write(17,*)
     write(17,'(a)') "unset xtics"
     write(17,'(a)') "unset ytics"
     write(17,'(a)') "set border"
!
! Put in functions for inequalities
!
     write(17,'(a11,f9.4,1x,a3,a3,f9.4,1x,a5,f9.4,a13)') "f1(x,y) = (", eqns(1,1),&
          &"* x", " + ", eqns(1,2), "* y >", eqns(1,nspecies_dep + 1), ")?f2(x,y):1/0"
     i=1
     if(nlimits > 0) then
        do i=2,nlimits+1
           if(i < 10 .and. i+1 < 10) then
              write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a5,f9.4,a3,i1,a9)') "f", i, "(x,y) = &
                   &(", eqns(i,1), "* x", " + ", eqns(i,2), "* y <", eqns(i,nspecies_dep &
                   &+ 1), ")?f", i+1, "(x,y):1/0"
           elseif(i < 10 .and. i+1 == 10) then
              write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a5,f9.4,a3,i2,a9)') "f", i, "(x,y) = &
                   &(", eqns(i,1), "* x", "+", eqns(i,2), "* y <", eqns(i,nspecies_dep &
                   &+ 1), ")?f", i+1, "(x,y):1/0"
           else
              write(17,'(a1,i2,a9,f9.4,a3,a3,f9.4,a5,f9.4,a3,i2,a9)') "f", i, "(x,y) = &
                   &(", eqns(i,1), "* x", "+", eqns(i,2), "* y <", eqns(i,nspecies_dep &
                   &+ 1), ")?f", i+1, "(x,y):1/0"
           endif
        enddo
     endif
     i = i-1
     if(i+1 < 10 .and. i+2 < 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a3,i1,a9)') "f", i+1, "(x,y) = (", eqns(i+1,1)&
             &, "* x", " < ", eqns(i+1,nspecies_dep + 1), ")?f", i+2, "(x,y):1/0"
     elseif(i+1 < 10 .and. i+2 == 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a3,i2,a9)') "f", i+1, "(x,y) = (", eqns(i+1,1)&
             &, "* x", " < ", eqns(i+1,nspecies_dep + 1), ")?f", i+2, "(x,y):1/0"
     else
        write(17,'(a1,i2,a9,f9.4,a3,a3,f9.4,a3,i2,a9)') "f", i+1, "(x,y) = (", eqns(i+1,1)&
             &, "* x", " < ", eqns(i+1,nspecies_dep + 1), ")?f", i+2, "(x,y):1/0"
     endif
     if(i+2 < 10 .and. i+3 < 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a3,i1,a9)') "f", i+2, "(x,y) = (", eqns(i+2,2)&
             &, "* y", " < ", eqns(i+2,nspecies_dep + 1), ")?f", i+3, "(x,y):1/0"
     elseif(i+2 < 10 .and. i+3 == 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a3,i2,a9)') "f", i+2, "(x,y) = (", eqns(i+2,2)&
             &, "* y", " < ", eqns(i+2,nspecies_dep + 1), ")?f", i+3, "(x,y):1/0"
     else
        write(17,'(a1,i2,a9,f9.4,a3,a3,f9.4,a3,i2,a9)') "f", i+2, "(x,y) = (", eqns(i+2,2)&
             &, "* y", " < ", eqns(i+2,nspecies_dep + 1), ")?f", i+3, "(x,y):1/0"
     endif
     if(i+3 < 10 .and. i+4 < 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a3,i1,a9)') "f", i+3, "(x,y) = (", eqns(i+3,1)&
             &, "* x", " > ", eqns(i+3,nspecies_dep + 1), ")?f", i+4, "(x,y):1/0"
     elseif(i+3 < 10 .and. i+4 == 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a3,i2,a9)') "f", i+3, "(x,y) = (", eqns(i+3,1)&
             &, "* x", " > ", eqns(i+3,nspecies_dep + 1), ")?f", i+4, "(x,y):1/0"
     else
        write(17,'(a1,i2,a9,f9.4,a3,a3,f9.4,a3,i2,a9)') "f", i+3, "(x,y) = (", eqns(i+3,1)&
             &, "* x", " > ", eqns(i+3,nspecies_dep + 1), ")?f", i+4, "(x,y):1/0"
     endif
     if(i+4 < 10) then
        write(17,'(a1,i1,a9,f9.4,a3,a3,f9.4,a7)') "f", i+4, "(x,y) = (", eqns(i+4,2), "* y"&
             &, " > ", eqns(i+4,nspecies_dep + 1), ")?1:1/0"
     else
        write(17,'(a1,i2,a9,f9.4,a3,a3,f9.4,a7)') "f", i+4, "(x,y) = (", eqns(i+4,2), "* y"&
             &, " > ", eqns(i+4,nspecies_dep + 1), ")?1:1/0"
     endif
     write(17,'(a)') "unset colorbox"
     write(17,'(a)') "set isosample 300, 300"
     write(17,'(a)') "set palette gray"
     write(17,'(a)') "unset xlabel"
     write(17,'(a)') "unset ylabel"
     write(17,'(a)') "set sample 500"
     write(17,'(a)') "set pm3d map"
     write(17,'(a7,f9.4,a5,f9.4,a19)') "splot [", eqns(1,nspecies_dep + 1), ":0] [", eqns(1,&
          &nspecies_dep + 1), ":0] f1(x,y) notitle"
     write(17,'(a)') "unset pm3d"
     write(17,'(a)') "set border 12"
     write(17,'(a)') "set key r c"
     write(17,'(a)') "set style fill transparent"
     write(17,'(a)') "set x2tics nomirror"
     write(17,'(a)') "set y2tics nomirror"
     write(17,'(a12,f9.4,a3)') "set xrange [", eqns(1,nspecies_dep + 1), ":0]"
     write(17,'(a12,f9.4,a3)') "set yrange [", eqns(1,nspecies_dep + 1), ":0]"
     write(17,'(a13,a2,a1)') "set x2label '", name(1), "'"
     write(17,'(a13,a2,a1)') "set y2label '", name(2), "'"
!
! Now write plot string - some competing phases may be left out if they have no definition
! in the 2D space. This can be the case when a large number of original elements is reduced
! by setting some chemical potentials to particular values. Effectively the coefficient can 
! become zero, and as some lines are defined by dividing by this coefficient, they must be 
! excluded
!
     lngstr = ""
!
! First write compound of interest formula
!
     write(lngstr,'(a4,f9.4,a3,f9.4,a8,a,a6)') "plot", eqns(1,nspecies_dep+1)/eqns(1,2), &
          &" - ", eqns(1,1)/eqns(1,2), " * x t '", trim(compounds(1)), "' lw 2"
!
! Now go through the competing phases and write for each one, skipping those that are undefined
!
     if(nlimits > 0) then
        do i=2,nlimits+1
           if(eqns(i,2) /= 0.d0) then
              write(string,'(a1,f9.4,a3,f9.4,a8,a,a6)') ",", eqns(i,nspecies_dep+1)/eqns(i,2), &
                   &" - ", eqns(i,1)/eqns(i,2), " * x t '", trim(compounds(i)), "' lw 2"
              lngstr = trim(lngstr)//trim(string)
           endif
        enddo
     endif
!
! Write plot string to file
!
     write(17,'(a)') trim(lngstr)
!
! Finish
!
     write(17,*)
     write(17,'(a)') "unset multiplot"
     close(unit=17)
     write(*,*)
     write(*,'(a)') "...done -> 2Dplot.plt can be loaded in gnuplot"
     write(*,'(a)') "        -> 2Dplot.dat can be pasted into Mathematica"
!
! Now write data to plot individual lines associated with the compounds to a .txt file, giving
! each compound separately, with two points defining the line. If a compound's line is parallel
! to the x or y axis, the line is defined according to the endpoints on the axes for the 
! compound of interest
!
     open(unit=17,status='replace',file="2Dplot.txt")
     write(17,'(a)') "#"//trim(compounds(1))
     write(17,'(2(f9.4,2x))') 0.d0, eqns(1,nspecies_dep+1)/eqns(1,2)
     write(17,'(2(f9.4,2x))') eqns(1,nspecies_dep+1)/eqns(1,1), 0.d0
     do i=2,nlimits+1
        write(17,*)
        write(17,'(a)') "#"//trim(compounds(i))
        if(eqns(i,1) == 0.d0) then
           write(17,'(2(f9.4,2x))') 0.d0, eqns(i,nspecies_dep+1)/eqns(i,2)
           write(17,'(2(f9.4,2x))') eqns(1,nspecies_dep+1)/eqns(1,1), &
                &eqns(i,nspecies_dep+1)/eqns(i,2)
        elseif(eqns(i,2) == 0.d0) then
           write(17,'(2(f9.4,2x))') eqns(i,nspecies_dep+1)/eqns(i,1), &
                &eqns(1,nspecies_dep+1)/eqns(1,2)
           write(17,'(2(f9.4,2x))') eqns(i,nspecies_dep+1)/eqns(i,1), 0.d0
        else
           if(eqns(i,nspecies_dep+1)/eqns(i,2) > 0.d0) then
              write(17,'(2(f9.4,2x))') eqns(1,nspecies_dep+1)/eqns(1,1), &
                   &eqns(i,nspecies_dep+1)/eqns(i,2) + -eqns(1,nspecies_dep+1)/eqns(1,1)&
                   &* eqns(i,1)/eqns(i,2)
           else
              write(17,'(2(f9.4,2x))') 0.d0, eqns(i,nspecies_dep+1)/eqns(i,2)
           endif
           if(eqns(i,nspecies_dep+1)/eqns(i,1) > 0.d0) then
              write(17,'(2(f9.4,2x))') (-eqns(1,nspecies_dep+1)/eqns(1,2) + &
                   &eqns(i,nspecies_dep+1)/eqns(i,2)) * eqns(i,2)/eqns(i,1), &
                   &eqns(1,nspecies_dep+1)/eqns(1,2)
           else
              write(17,'(2(f9.4,2x))') eqns(i,nspecies_dep+1)/eqns(i,1), 0.d0
           endif
        endif
     enddo
!
! Finish
!
     close(unit=17)
     write(*,*)
     write(*,'(a)') "        -> 2Dplot.txt contains points to plot lines"
!
  else
     continue
  endif
!
! Tell user how to do a restart and fix one of the chemical potentials
!
  if(nspecies >= 3) then
     write(*,*)
     write(*,'(a)') "Copy 'results.dat' to 'restart.dat' in order to run the"
     write(*,'(a)') "calculation again with a different dependent variable, with"
     write(*,'(a)') "additional competing phases, or with one of the chemical"
     write(*,'(a)') "potentials set at a fixed value."
     write(*,*)
  else
     write(*,*)
     write(*,'(a)') "Copy 'results.dat' to 'restart.dat' in order to run the"
     write(*,'(a)') "calculation again with a different dependent variable, or"
     write(*,'(a)') "with additional competing phases."
     write(*,*)
  endif
!
200 continue
!
  close(unit=11)
  close(unit=12)
!
end program cplap
!
!
!
!
subroutine all_eqns_solve(eqns,dim1,dim2,results,nlim,numres)
!
! This subroutine reads in an array 'eqns' containing a number dim1
! of linear equations, with dim2-1 unknowns in each equation. The 
! routine forms all possible combinations of the linear equations
! and solves each system (if a solution exists), using LU decomposition
! and back substitution (using modified versions of the numerical
! recipes routines ludcmp and lubksb (Numerical Recipes in Fortran 77. 
! The Art of Scientific Computing, 2nd Edition, 1992 W. H. Press et al.)
! The solutions are returned in the array 'results'
!
! J. Buckeridge October 2011
!
! HACK - 10/03/2016: Added check to see if solution found is compatible with 
!                    limiting inequalities. Previously this check was done
!                    outside the subroutine, but putting it here saves space.
!
  implicit none
!
  integer, parameter :: nmax = 300
  real*8, parameter :: small = 1.d-12
!
  integer i, j, k, dim1, dim2, num, indx(nmax), ival(nmax)
  integer, intent(inout) :: numres
  integer, intent(inout) :: nlim
!  real*8, intent(inout) :: eqns(:,:)
  !real*8, dimension(:,:), intent(inout) :: results
  !real*8 eqns(dim1,dim2), results(nmax*nmax,dim2)
  real*8 eqns(dim1,dim2), sum1, results(nmax*nmax,dim2)
  real*8 system(size(eqns,2),size(eqns,2)), vector(size(eqns,2)), d
  logical check1, checksol1, checksol2
  character*20 string
!
! Set num equal to the dimensions of the linear system to be solved (i.e. 
! equal to dim2 - 1)
!
  num = dim2 - 1
!
! Initialize numres as zero
!
  numres = 0
!
! There will be dim1 C num systems to be solved (if solutions for each
! particular system exist). Go through each system and attempt to solve.
! This is done using an array of length num, each entry of which is
! increased in the correct order so that all combinations of length num
! from the array of length dim1 are taken
!
  do i=1,num
     ival(i) = i
  enddo
100 continue
!
! Set the num x num array to be solved, and the vector in which the
! solutions are returned from the numerical recipes routines. Re-zero
! both first
!
  system = 0.d0
  vector = 0.d0
  do i=1,num
     do j=1,num
        system(i,j) = eqns(ival(i),j)
     enddo
     vector(i) = eqns(ival(i),dim2)
  enddo
!
! Call numerical recipes routine to get LU decomposition. If this fails
! then there is no solution for this system of equations.
!
! If solution is found, check its compatibility with the limiting inequalities
! Write solutions to scratch file to save dynamical memory
!
  check1 = .true.
  call ludcmp(system(1:num,1:num),indx(1:num),d,num,check1)
  if(check1) then
     call lubksb(system(1:num,1:num),indx(1:num),vector(1:num),num)
     if(any(vector >= (1.d0/small))) then
        goto 102
     else
        checksol1 = .true.
        do j=1,dim1
           sum1 = 0.d0
           do k=1,num
              sum1 = sum1 + vector(k) * eqns(j,k)
           enddo
           if(j == 1 .or. j > nlim + num + 1) then
              if(sum1 - eqns(j,dim2) < -small) then
                 checksol1 = .false.
              endif
           else
              if(sum1 - eqns(j,dim2) > small) then
                 checksol1 = .false.
              endif
           endif
        enddo
        if(checksol1) then
           numres = numres + 1
           if(numres > nmax*nmax) then
              write(*,'(a)') "*** sorry, too much memory used ***"
              stop
           endif
           results(numres,1:num) = vector(1:num)
        endif
     endif
  endif
102 continue
!
  i = 0
  do
     ival(num-i) = ival(num-i)+1
     if(ival(num-i) <= dim1-i) then
        goto 100
     else
        i = i+1
        if(i >= num) goto 101
        do j=1,i
           if(j == 1) then
              ival(num-i+j) = ival(num-i+j-1)+2
           else
              ival(num-i+j) = ival(num-i+j-1)+1
           endif
        enddo
     endif
  enddo
101 continue
!
end subroutine all_eqns_solve
!
!
!
!
SUBROUTINE ludcmp(a,indx,d,nn,check)
!
! HACK: replace all references to modules with own code. This includes replacing
! function calls with code blocks for the actual function. Done to avoid SP/DP etc
! confusion
! HACK: change precision of variables to maximum
! HACK: input nn to declare array sizes
! HACK: change declarations to F77 style
!
!
! HACK (oct 2011) : Added check, to check if matrix is singular (if so no solution
!                   will be found)
!
!
 !USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
 IMPLICIT NONE
 integer nn, indx(nn), j, n, imax
 real*8 a(nn,nn), d, vv(nn)
 logical, intent(inout) :: check
 !REAL*8, DIMENSION(:,:), INTENT(INOUT) :: a
 !REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
 !integer indx(2*nn)
 !INTEGER, DIMENSION(:), INTENT(OUT) :: indx
 !INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
 !REAL*8, INTENT(OUT) :: d
 !REAL(SP), INTENT(OUT) :: d
!Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!output vector of length N that records the row permutation effected by the partial pivoting;
!d is output as ±1 depending on whether the number of row interchanges was even or odd,
!respectively. This routine is used in combination with lubksb to solve linear equations or
!invert a matrix.
 !real*8 vv(2*nn)
 !REAL*8, DIMENSION(size(a,1)) :: vv
 !REAL(SP), DIMENSION(size(a,1)) :: vv
!vv stores the implicit scaling of each row.
 REAL*8, PARAMETER :: TINY=1.0d-20
 !REAL(SP), PARAMETER :: TINY=1.0e-20_sp
!A small number.
 !INTEGER :: j,n,imax
 !INTEGER(I4B) :: j,n,imax
!
! HACK: dummy array for imaxloc function
!
 integer, dimension(1) :: imaxloc_r
!
! HACK: dummy array for swap function
!
 real*8, dimension(size(a,1)) :: dum
!
! HACK: dummy array for outerprod function
!
 REAL*8, DIMENSION(size(a),size(a)) :: outerprod_r
!
! HACK: replace assert_eq with code block
!
!n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
 if (size(a,1) == size(a,2) .and. size(a,2) == size(indx)) then
    n=size(a,1)
 else
    write (*,*) 'nrerror: an assert_eq failed with this tag: ludcmp'
    STOP 'program terminated by assert_eq3'
 end if
 d=1.d0
 !d=1.0
!No row interchanges yet.
 vv=maxval(abs(a),dim=2)
!Loop over rows to get the implicit scaling
!
! HACK: replace call to nrerror
!
!
! ****HACK (oct 2011) : if matrix is singular return check=true
!                       (also change test to vv <= TINY)
!
! if (any(vv == 0.0)) then
 if (any(vv <= TINY)) then
    check = .false.
    return
 endif
!write(*,*) "nrerror: singular matrix in ludcmp"
 !if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
!information.
!There is a row of zeros.
!
! HACK: change to double precision
!
 vv=1.d0/vv
 !vv=1.0_sp/vv
!Save the scaling.
 do j=1,n
!
! HACK: replace imacloc with code block
!
    imaxloc_r = maxloc(vv(j:n)*abs(a(j:n,j)))
    imax=(j-1)+imaxloc_r(1)
    !imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
!Find the pivot row.
    if (j /= imax) then
!Do we need to interchange rows?
!Yes, do so...
!
! HACK: replace swap with code block
!
       dum=a(imax,:)
       a(imax,:)=a(j,:)
       a(j,:)=dum
       d=-d
!call swap(a(imax,:),a(j,:))
!...and change the parity of d.
       vv(imax)=vv(j)
!Also interchange the scale factor.
    end if
    indx(j)=imax
    if (a(j,j) == 0.0) a(j,j)=TINY
!If the pivot element is zero the matrix is singular (at least to the precision of the al-
!gorithm). For some applications on singular matrices, it is desirable to substitute TINY
!for zero.
    a(j+1:n,j)=a(j+1:n,j)/a(j,j)
!Divide by the pivot element.
!
! HACK: replace outerprod with code block
!
    outerprod_r = 0.d0
    outerprod_r(1:size(a(j+1:n,j)),1:size(a(j,j+1:n))) = spread(a(j+1:n,j),&
       &dim=2,ncopies=size(a(j,j+1:n))) * spread(a(j,j+1:n),dim=1,ncopies=size(a(j+1:n,j)))
    a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod_r(1:size(a(j+1:n,j)),1:size(a(j,j+1:n)))
    !a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!Reduce remaining submatrix.
 end do
END SUBROUTINE ludcmp
!
!
!
!
SUBROUTINE lubksb(a,indx,b,nn)
!
! HACK: replace all references to modules with own code. This includes replacing
! function calls with code blocks for the actual function. Done to avoid SP/DP etc
! confusion
! HACK: change precision of variables to maximum
! HACK: input nn to declare array sizes
! HACK: change declarations to F77 style
!
 !USE nrtype; USE nrutil, ONLY : assert_eq
 IMPLICIT NONE
 integer nn, indx(nn), i, n, ii, ll
 real*8 a(nn,nn), b(nn), summ
 !REAL*8, DIMENSION(:,:), INTENT(IN) :: a
 !REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
 !integer indx(2*nn)
 !INTEGER, DIMENSION(:), INTENT(IN) :: indx
 !INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
 !REAL, DIMENSION(:), INTENT(INOUT) :: b
 !REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
!Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!as the original matrix A, but rather as its LU decomposition, determined by the routine
!ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!input as the right-hand-side vector B, also of length N , and returns with the solution vector
!X. a and indx are not modified by this routine and can be left in place for successive calls
!with different right-hand sides b. This routine takes into account the possibility that b will
!begin with many zero elements, so it is efficient for use in matrix inversion.
 !INTEGER :: i,n,ii,ll
 !INTEGER(I4B) :: i,n,ii,ll
 !REAL*8 :: summ
 !REAL(SP) :: summ
!
! HACK: replace assert_eq with code block
!
!n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
 if (size(a,1) == size(a,2) .and. size(a,2) == size(indx)) then
    n=size(a,1)
 else
    write (*,*) 'nrerror: an assert_eq failed with this tag: lubksb'
    STOP 'program terminated by assert_eq3'
 end if
 ii=0
!When ii is set to a positive value, it will become the in-
!dex of the first nonvanishing element of b. We now do
 do i=1,n
!the forward substitution, equation (2.3.6). The only new
    ll=indx(i)
!wrinkle is to unscramble the permutation as we go.
    summ=b(ll)
    b(ll)=b(i)
    if (ii /= 0) then
       summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
    else if (summ /= 0.0) then
       ii=i
!A nonzero element was encountered, so from now on we will
    end if
!have to do the dot product above.
    b(i)=summ
 end do
 do i=n,1,-1
!Now we do the backsubstitution, equation (2.3.7).
    b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
 end do
END SUBROUTINE lubksb
!
