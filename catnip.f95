      program catnip

      !use IFPORT
! comment ^ out if not using ifort
! ifort -Tf catnip.f95 -free -Ofast -o catnip.e
! gfortran catnip.f95 -O3 -o catnip.e
      implicit none

      logical:: help,ential
      real, allocatable:: xyz(:,:),nxyz(:,:),tryxyz(:,:)
      real:: dist,c1,c2,com1,com2,box(3),bot,top,cx,cy,cz,randumb,psi,&
      a,b,c,vx,vy,vz,theta,dxyz(3),a1,a2,a3,b1,b2,b3,nb,pie,whattimeisit
      double precision:: R(3,3),dx,dy,dz
      integer:: i,j,k,l,m,n,o,tries,attempt,boxatoms,newatoms,lim,&
      ilist(8),deleted,groups,molsofthat,nextatom,nextmol,igroup,&
      rngesus,ios,intnum,moldel,lugicul
      integer, allocatable:: groupies(:),molnum(:),atnum(:),nmolnum(:)& 
      ,natnum(:),catdel(:),tempint(:)
      character:: str*99,agroup*90,g1*4,g2*4,outbox*90,outtop*90,&
      molfile*90,boxfile*90,dlist*90,ndxfile*90,axis,topfile*90,fi*3&
      ,fo*3,fn*3
      character, allocatable:: molres(:)*4,atres(:)*4,nmolres(:)*4,&
      natres(:)*4,delgroups(:)*4

      help=.False.

      if(iargc().eq.0) help=.True.
! help if there is no info at runtime

      i=0
      igroup=-1
      boxfile="a1"
      topfile="a1"
      outbox="catnip.gro"
      outtop="catnip.top"
! default outs
      dlist=" "
      ndxfile="index.ndx"
      agroup="-1"
      tries=-1
      g1="6.02"
      g2="3.14"
      axis='z'

      do i=1,iargc()
        call getarg(i,str)
        if(str.eq."-c") then
          call getarg(i+1,boxfile)
        endif
        if(str.eq."-oc") then
         call getarg(i+1,outbox)
        endif
        if(str.eq."-p") then
         call getarg(i+1,topfile)
        endif
        if(str.eq."-op") then
          call getarg(i+1,outtop)
        endif
        if(str.eq."-ci") then
          call getarg(i+1,molfile)
        endif
        if(str.eq."-dl") then
          call getarg(i+1,dlist)
        endif
        if(str.eq."-n") then
         call getarg(i+1,str)
         read(str,'(i90)') tries
        endif
        if(str.eq."-g") then
          call getarg(i+1,agroup)
          if(agroup.eq."interface") then
            call getarg(i+2,g1)
            call getarg(i+3,g2)
            if(i+4.le.iargc()) call getarg(i+4,axis)
          endif
        endif
        if(str.eq."-ndx") then
          call getarg(i+1,ndxfile)
        endif
        if(str.eq."-h") then
          help=.True.
          exit
        endif
      enddo
! read command line stuff 

      open(1,file="/dev/stdout")

      if(help) then
        str=adjustl("")
        write(1,9) str
        str=adjustl("       catnip is used to put residues from one &
          &file into another")
        write(1,9) str
        str=adjustl("")
        write(1,9) str
        str=adjustl("-c   : input system coordinate file (gro/pdb) to &
          &which molecules are added")
        write(1,9) str
        str=adjustl("-oc  : [ optional | catnip.gro ] output system &
          &coordinate file (gro/pdb)")
        write(1,9) str
        str=adjustl("-p   : input topology file (.top) which &
          &corresponds to the system coordinate file")
        write(1,9) str
        str=adjustl("-op  : [ optional | catnip.top ] output topology &
          &file (.top)")
        write(1,9) str
        str=adjustl("-ci  : input molecule file (gro/pdb) which is put &
          &in the system file")
        write(1,9) str
        str=adjustl("-dl  : [ optional ] comma-separated list or file &
          &containing removable residues (e.g. LIG,DRG)")
        write(1,9) str
        str=adjustl("-n   : [ optional | prompted ] number of new &
          &molecules to put into the system")
        write(1,9) str
        str=adjustl("-g   : [ optional | prompted ] group name, number,&
          & or 'interface RES1 (or ceil) RES2 x||y||z'")
        write(1,9) str
        str=adjustl("-ndx : [ optional | index.ndx ] specify index &
          &file")
        write(1,9) str
        str=adjustl("-h   : help")
        write(1,9) str
        goto 97
      endif 
! useful info ^
      
      pie=acos(-1.)
 
      do while(tries.lt.0) 
        write(1,*) "How many to add?"
        read(5,*) tries
      enddo
! if not given a command line option, it will ask how many to add

      inquire(file=dlist,exist=ential)

      if(ential) then

        n=0
        open(12,file=dlist)
        do
          read(12,*,iostat=ios) str
          if(ios.ne.0) exit
          n=n+1
        enddo

        allocate(delgroups(n))

        rewind(12)

        do i=1,n
          read(12,*) delgroups(i)
        enddo

        close(12)
! read file for removable molecules
      else

        n=1
        do k=1,90
          if(dlist(k:k).eq.",") then
            n=n+1
          endif
        enddo
! split list on commas  
        allocate(delgroups(n))

        n=1
        l=1
        do i=1,90
          if(dlist(i:i).eq.",") then
            delgroups(n)=dlist(l:i-1)
            n=n+1
            l=i+1
          endif
        enddo

        delgroups(n)=dlist(l:90)
      
      endif
! read commandline delete stuff stuff

      write(1,*) "removable molecules:",delgroups

      inquire(file=ndxfile,exist=ential)

      if(.not.ential) then
        write(1,*)"generating index file"
        str="echo q | gmx -quiet make_ndx -f "//trim(boxfile)//"&
          &&>/dev/null"
        lugicul=system(str)
        ndxfile="index.ndx"
      else
        write(1,*) "using pre-existing index file"
      endif
! make new index if there isn't one

      inquire(file=boxfile,exist=ential)
      if(.not.ential) then
       write(1,*) "box file not properly specified, use -h for help"
       goto 97
      endif
! goto 97 ends the program

      boxfile=adjustr(boxfile)
      fi=boxfile(88:90)
      boxfile=adjustl(boxfile)

      !inquire(file=topfile,exist=ential)
      !if(.not.ential) then
      ! write(1,*) "topology file not properly specified, &
      ! &use -h for help"
      ! goto 97
      !endif

      write(1,*) "reading index file..."

      open(6,file=ndxfile)

      groups=0

37    if(verify(trim(agroup),'0123456789').ne.0.and.agroup.ne."-1") then
        agroup='[ '//trim(agroup)//' ]'
      else
        read(agroup,'(i90)') igroup
      endif
! check if the user gave a number or a name for the group, and search
! according to their input
      l=0
      do
        read(6,'(a)',iostat=ios) str
        if(ios.ne.0) exit
        str=adjustl(str)
        if(str(1:1).eq.'[') then
          if(groups.eq.igroup.or.str.eq.agroup) then
            do
              read(6,'(a)',iostat=ios) str
              str=adjustl(str)
              if(ios.ne.0.or.str(1:1).eq.'[') exit
              n=90
              str=' '//str
              do k=89,1,-1
                if(str(k:k).eq.' '.and.str(k+1:k+1).ne.' ') then
                  read(str(k:n),'(i6)') m
                  n=k
                  if(allocated(groupies)) then
                    intnum=size(groupies)
                    allocate(tempint(intnum+1))
                    do o=1,intnum
                      tempint(o)=groupies(o)
                    enddo
                    tempint(intnum+1)=m
                    deallocate(groupies)
                    call move_alloc(tempint,groupies)
                  else
                    allocate(groupies(1))
                    groupies(1)=m
                  endif
! append each number to the list
! random item in this list will be chosen later
                  j=k
                endif
              enddo
            enddo
          else
            if(agroup.eq.'-1'.and.igroup.eq.-1) then
              write(1,*) groups,trim(str)
! write groups if the group hasn't been changed from the initial
            endif
            groups=groups+1
! count the groups
          endif
        endif
      enddo
! read index file and members of the group specified

      if(.not.allocated(groupies).and.agroup.ne.'-1') then
        write(1,*) trim(agroup)," not found in index"
      endif
! say if the selected group wasn't found
      if(agroup.eq."-1") then
        write(1,*) "Which group?"
! ask if group is not found in index
        read(*,*) agroup
        if(agroup.eq."interface") then
          write(1,*) "interface group 1:"
          read(*,*) g1
          write(1,*) "interface group 2:"
          read(*,*) g2
          write(1,*) "axis:"
          read(*,*) axis
          agroup="[ interface ]"
        else
          groups=0
          rewind(6)
          goto 37
! go back and read index
        endif
      endif
! ask about group is not specified or not in index

      close(6)

      open(7,file=boxfile)
      
      if(fi.eq.'gro') then
        read(7,'(a)') str
        read(7,*) boxatoms
      else
        boxatoms=0
        do
          read(7,'(a)',iostat=ios) str
          if(ios.ne.0) exit
          if(str(1:4).eq.'ATOM'.or.str(1:6).eq.'HETATM') then
            boxatoms=boxatoms+1
          endif
        enddo
        rewind(7)
      endif
      
      inquire(file=molfile,exist=ential)
      if(.not.ential) then
        write(1,*) "input file not properly specified, use -h for help"
        goto 97
      endif

      molfile=adjustr(molfile)
      fn=molfile(88:90)
      molfile=adjustl(molfile)

      open(8,file=molfile)
      
      if(fn.eq.'gro') then
        read(8,*) str
        read(8,*) newatoms
      else
        newatoms=0
        do
          read(8,'(a)',iostat=ios) str
          if(ios.ne.0) exit
          if(str(1:4).eq.'ATOM'.or.str(1:6).eq.'HETATM') then
            newatoms=newatoms+1
          endif
        enddo
        rewind(8)
      endif

      allocate(molnum(boxatoms+tries*newatoms),molres(boxatoms+tries*&
       newatoms),atres(boxatoms+tries*newatoms),atnum(boxatoms+tries*&
       newatoms),xyz(boxatoms+tries*newatoms,3),catdel(boxatoms))

      catdel(:)=-1
! catdel is the list of deleted molecules
      if(fi.eq.'gro') then
        do i=1,boxatoms
          read(7,99) molnum(i),molres(i),atres(i),atnum(i),&
            (xyz(i,j),j=1,3)
! read box file
        enddo
      read(7,*) box
      else
        box=(/0.,0.,0./)
        i=0
        do
          read(7,'(a)',iostat=ios) str
          if(ios.ne.0) exit
          if(str(1:4).eq.'ATOM'.or.str(1:6).eq.'HETATM') then
            i=i+1
            read(str,96) atnum(i),atres(i),molres(i),molnum(i),&
              (xyz(i,j),j=1,3)
            do j=1,3
              xyz(i,j)=xyz(i,j)/10.
            enddo
          endif
        if(str(1:6).eq.'CRYST1') then
          backspace(7)
          read(7,*) str,box
          box=box/10.
        endif
        enddo
        if(any(box.eq.0.)) then
          box(1)=maxval(xyz(:,1))-minval(xyz(:,1))
          box(2)=maxval(xyz(:,2))-minval(xyz(:,2))
          box(3)=maxval(xyz(:,3))-minval(xyz(:,3))
        endif
      endif
      close(7)

      if(verify(axis,'xyz').ne.0) then
        write(1,*) "enter axis (x, y, or z)"
        read(*,*) axis
! ask about the axis if it isn't x, y, or z
      endif

      if(agroup.eq."[ interface ]") then
        write(1,*) "determining interface..."
        if(axis.eq."x") lim=1
        if(axis.eq."y") lim=2
        if(axis.eq."z") lim=3
    
        c1=0.
        c2=0.
        com1=0.
        com2=0.
        do i=1,boxatoms
          if(molres(i).eq.g1) then
            com1=com1+xyz(i,lim)
            c1=c1+1.
          endif
          if(molres(i).eq.g2) then
            com2=com2+xyz(i,lim)
            c2=c2+1.
          endif
        enddo
     
        if(c1.gt.0.) then
          com1=com1/c1
        else
          com1=0.
        endif
  
        if(c2.gt.0.) then
          com2=com2/c2
        else
          com2=0.
        endif
  
        if(g1.eq.'ceil') com1=box(lim)
        if(g2.eq.'ceil') com2=box(lim)
  
        top=max(com1,com2)
        bot=min(com1,com2)
  
        do i=1,boxatoms
          if(xyz(i,lim).lt.top.and.xyz(i,lim).gt.bot) then
            if(allocated(groupies)) then
              intnum=size(groupies)
              allocate(tempint(intnum+1))
              do o=1,intnum
                tempint(o)=groupies(o)
              enddo
              tempint(intnum+1)=atnum(i)
              deallocate(groupies)
              call move_alloc(tempint,groupies)
            else
              allocate(groupies(1))
              groupies(1)=atnum(i)
           endif
          endif
        enddo
      endif
! find interface if needed
! 'interface' is the region between the middles of the selected groups
! ceil is the top of the box
  
      allocate(nmolnum(newatoms),nmolres(newatoms),natres(newatoms),& 
        natnum(newatoms),nxyz(newatoms,3),tryxyz(newatoms,3))

      cx=0.
      cy=0.
      cz=0.

      if(fn.eq.'gro') then
        do l=1,newatoms
          read(8,99) nmolnum(l),nmolres(l),natres(l),natnum(l), &
            (nxyz(l,j),j=1,3)
          cx=cx+nxyz(l,1)/real(newatoms)
          cy=cy+nxyz(l,2)/real(newatoms)
          cz=cz+nxyz(l,3)/real(newatoms)
! read molecules to add and find the center
        enddo
      else
        l=0
        do
          read(8,'(a)',iostat=ios) str
          if(ios.ne.0) exit
          if(str(1:4).eq.'ATOM'.or.str(1:6).eq.'HETATM') then
            l=l+1
            read(str,96) natnum(l),natres(l),nmolres(l),nmolnum(l),&
              (nxyz(l,j),j=1,3)
              do j=1,3
                nxyz(l,j)=nxyz(l,j)/10.
              enddo
            cx=cx+nxyz(l,1)/real(newatoms)
            cy=cy+nxyz(l,2)/real(newatoms)
            cz=cz+nxyz(l,3)/real(newatoms)
          endif
        enddo
      endif

      close(8)

      do l=1,newatoms
        nxyz(l,1)=nxyz(l,1)-cx
        nxyz(l,2)=nxyz(l,2)-cy
        nxyz(l,3)=nxyz(l,3)-cz
      enddo
! move new molecules relative to the center(ish)
      
      moldel=0
      deleted=0

      do k=1,tries

        call date_and_time(VALUES=ilist)

        call srand(ilist(8)+ilist(7)+deleted+k)
! use current time to be more randumb

        randumb=rand()
        randumb=rand()
! the 1st randumb number is randumb, always low

        rngesus=groupies(floor((size(groupies)*randumb)+1.))
! pick the chosen one

        write(1,*) "placing new molecule",k,"by",rngesus

        do l=1,newatoms
          tryxyz(l,:)=nxyz(l,:)+xyz(rngesus,:)
        enddo
! place new molecules on top of existing ones 

        attempt=0

        write(1,*) "adjusting orientation..."
        
        i=0
        outer: do
          i=i+1
          if(i.gt.boxatoms) exit
          if(.not.any(delgroups.eq.molres(i))) then
            l=0
            do 
              l=l+1
              if(l.gt.newatoms) exit
              call sqdist(xyz(i,:),tryxyz(l,:),box,dist) 
              if(dist.lt.0.05) then
                i=0
                l=0
                attempt=attempt+1
! if the molecules that are being added are too close to something that
! isn't allowed to be removed, it will be translated or rotated to try
! to get them farther apart
                call sqdist(xyz(i,:),xyz(rngesus,:),box,a)
                call sqdist(tryxyz(l,:),xyz(rngesus,:),box,b)
                call sqdist(xyz(i,:),tryxyz(l,:),box,c)
  
                if(a.gt.1..and.b.gt.1.) then
                  psi=acos((c-a-b)/(-2.*sqrt(a*b)))
                else
                  psi=0.
                endif
 
                if(abs(psi).gt.1.2.or.psi.eq.0..or.c.lt..01) then
                  if(abs(xyz(i,1)-tryxyz(l,1)).gt.abs(box(1)-&
                    abs(xyz(i,1)-tryxyz(l,1)))) then
                    if(xyz(i,1).gt.tryxyz(l,1)) then
                      dxyz(1)=box(1)-abs(xyz(i,1)-tryxyz(l,1))
                    else
                      dxyz(1)=(tryxyz(l,1)-xyz(i,1))-box(1)
                    endif
                  else
                    dxyz(1)=tryxyz(l,1)-xyz(i,1)
                  endif
                  if(abs(xyz(i,2)-tryxyz(l,2)).gt.abs(box(2)-&
                    abs(xyz(i,2)-tryxyz(l,2)))) then
                    if(xyz(i,2).gt.tryxyz(l,2)) then
                      dxyz(2)=box(2)-abs(xyz(i,2)-tryxyz(l,2))
                    else
                      dxyz(2)=(tryxyz(l,2)-xyz(i,2))-box(2)
                    endif
                  else
                    dxyz(2)=tryxyz(l,2)-xyz(i,2)
                  endif
                  if(abs(xyz(i,3)-tryxyz(l,3)).gt.abs(box(3)-&
                    abs(xyz(i,3)-tryxyz(l,3)))) then
                    if(xyz(i,3).gt.tryxyz(l,3)) then
                      dxyz(3)=box(3)-abs(xyz(i,3)-tryxyz(l,3))
                    else
                      dxyz(3)=(tryxyz(l,3)-xyz(i,3))-box(3)
                    endif
                  else
                    dxyz(3)=tryxyz(l,3)-xyz(i,3)
                  endif
                  nb=0
                  do j=1,3
                    nb=nb+dxyz(j)**2
                  enddo
                  if(n.ne.0.) then
                    do j=1,3
                     dxyz(j)=dxyz(j)/(4.*nb)
                    enddo
                  else
                    dxyz=(/0.25,0.0,0.0/)
                  endif
          
                  do j=1,newatoms
                    tryxyz(j,:)=tryxyz(j,:)+dxyz
! translate it away
                  enddo
                else
                  a1=xyz(i,1)-tryxyz(l,1)
                  a2=xyz(i,2)-tryxyz(l,2)
                  a3=xyz(i,3)-tryxyz(l,3)
                  b1=xyz(i,1)-xyz(rngesus,1)
                  b2=xyz(i,2)-xyz(rngesus,2)
                  b3=xyz(i,3)-xyz(rngesus,3)
        
                  vx=(a2*b3)-(a3*b2)
                  vy=-(a1*b3)+(a3*b1)
                  vz=(a1*b2)-(a2*b1)
                  nb=sqrt(vx**2+vy**2+vz**2)
                  if(nb.ne.0.) then
                    vx=vx/nb
                    vy=vy/nb
                    vz=vz/nb
                  else
                    vx=1.0
                    vy=0.0
                    vz=0.0
                  endif
! find unit vector parallel to axis of rotation
! if the atoms are on top of each other, use the x axis

                  theta=0.5
! rotation angle is hardcoded for now

                  R(1,1)=cos(theta)+vx*vx*(1-cos(theta))
                  R(1,2)=vx*vy*(1-cos(theta))-vz*sin(theta)
                  R(1,3)=vx*vz*(1-cos(theta))+vy*sin(theta)
                  R(2,1)=vx*vy*(1-cos(theta))+vz*sin(theta)
                  R(2,2)=cos(theta)+vy*vy*(1-cos(theta))
                  R(2,3)=vy*vz*(1-cos(theta))-vx*sin(theta)
                  R(3,1)=vx*vz*(1-cos(theta))-vy*sin(theta)
                  R(3,2)=vy*vz*(1-cos(theta))+vx*sin(theta)
                  R(3,3)=cos(theta)+vz*vz*(1-cos(theta))
! rotation matrix

                  do j=1,newatoms
                    dx=nxyz(j,1)*R(1,1)+nxyz(j,2)*R(1,2)+nxyz(j,3)*&
                      R(1,3)
                    dy=nxyz(j,1)*R(2,1)+nxyz(j,2)*R(2,2)+nxyz(j,3)*&
                      R(2,3)
                    dz=nxyz(j,1)*R(3,1)+nxyz(j,2)*R(3,2)+nxyz(j,3)*&
                      R(3,3)
                    nxyz(j,:)=(/dx,dy,dz/)
                  enddo    

                  do j=1,newatoms
                    tryxyz(j,:)=nxyz(j,:)+xyz(rngesus,:)
! rotate it
                  enddo
                endif
              endif
              if(attempt.ge.30) exit outer
            enddo
          endif
        enddo outer
! move molecules away from ones that are too close, but can't be removed
! it's got 30 tries (rotations and/or translations) to move the molecule away

        write(1,*) "removing molecules..."

        do i=1,boxatoms
          if(any(delgroups.eq.molres(i)).and.&
            .not.any(catdel.eq.molnum(i))) then
! user-specified molecules can be removed if they are too close to the
! new molecules
            do l=1,newatoms
              call sqdist(xyz(i,:),tryxyz(l,:),box,dist)
              if(dist.lt.0.02) then
                moldel=moldel+1
                catdel(moldel)=molnum(i)
                deleted=deleted+count(molnum(i).eq.molnum(1:boxatoms))
! keep track of number of atoms removed and which molecules will be
! removed
                exit
              endif
            enddo
          endif
        enddo
! delete molecules that are too close and can be removed

        do l=1,newatoms
          molnum(boxatoms+l)=molnum(boxatoms)+nmolnum(l)
          molres(boxatoms+l)=nmolres(l)
          atres(boxatoms+l)=natres(l)
          xyz(boxatoms+l,:)=tryxyz(l,:)
        enddo
! add the new stuff to the list
     
        boxatoms=boxatoms+newatoms         
! update number of atoms in the box (ignoring removed ones for now)
      
      enddo
  
      write(1,*) "writing output files..."

      open(9,file=outbox)
      outbox=adjustr(outbox)
      fo=outbox(88:90)
      outbox=adjustl(outbox)
      if(fo.eq.'pdb') then
        box=box*10.
        xyz=xyz*10.
        atres=adjustl(atres)
        molres=adjustl(molres)
      endif
      open(11,file=outtop)
      
      inquire(file=topfile,exist=ential)
      if(ential) then
        open(10,file=topfile)

        do
          read(10,'(a)',iostat=ios) str
          write(11,'(a)') str
          str=adjustl(str)
          if(str.eq."[ molecules ]".or.ios.ne.0) exit
        enddo

! copy everything from the old topology up until [ molecules ]

        close(10)
      else
        write(11,'(a)') '[ molecules ]'
      endif
      
      if(fo.eq.'gro') then
        write(9,'(a)') "catnip conffile"
        write(9,'(i5)') boxatoms-deleted
      endif

      if(any(catdel.eq.molnum(1))) then
        nextmol=0
      else
        nextmol=1
      endif
      nextatom=0
      molsofthat=0

      do i=1,boxatoms
        if(.not.any(catdel.eq.molnum(i))) then
          nextatom=nextatom+1
          if(i.gt.1) then 
            if(molnum(i).ne.molnum(i-1)) then
              nextmol=nextmol+1
              molsofthat=molsofthat+1
              if(fo.eq.'pdb') write(9,'(a3,3x,i5,6x,a4)') &
                'TER',nextatom,molres(i)
            endif
! renumber atoms and molecules so they are consecutive
          endif
! count the new number of things new number of things to topology
! write new box file
          if(fo.eq.'gro') then
            write(9,99) & 
              nextmol,molres(i),atres(i),nextatom,(xyz(i,j),j=1,3)
          else
            write(9,98) "HETATM",nextatom,atres(i),molres(i),nextmol&
            ,(xyz(i,j),j=1,3),1.00,0.00 
          endif
        endif
        if(molres(i).ne.molres(i-1).and.i.gt.1) then
          write(11,'(a4,2x,i5)') molres(i-1),molsofthat
          molsofthat=0
        endif
      enddo
      if(.not.any(catdel.eq.molnum(boxatoms)).and.nextmol.gt.1) &
        molsofthat=molsofthat+1
      write(11,'(a4,2x,i5)') molres(size(molres)),molsofthat
      if(fo.eq.'gro') write(9,'(3(2x,f8.5))') box
      if(fo.eq.'pdb') write(9,'(a)') "END"
! write all atoms to the new files

      close(9)
      close(11)

97    call cpu_time(whattimeisit)
! 97 is the go-to for exiting the program early for errors and stuff

      if(allocated(groupies)) deallocate(groupies)
      if(allocated(molnum)) deallocate(molnum,molres,atres,atnum,xyz,&
        &catdel)
      if(allocated(nmolnum)) deallocate(nmolnum,nmolres,natres,natnum,&
        nxyz,tryxyz)
      if(allocated(delgroups)) deallocate(delgroups)
! make sure things are deallocated, even if things are messed up

      write(1,'(a4,1x,f8.5,a1)') 'took',whattimeisit,'s'

9     format(2x,a99)
96    format(6x,i5,2x,a4,a4,1x,i4,4x,3(1x,f7.3))
98    format(a6,i5,2x,a4,a4,1x,i4,4x,3(1x,f7.3),2x,f4.2,1x,f4.2)
99    format(i5,a4,2x,a4,i5,3(2x,f6.3))

      close(1)

      end program catnip

      subroutine sqdist(xyz1,xyz2,box,dist2)
! squared distance between atoms
      implicit none

      real, intent(in):: xyz1(3),xyz2(3),box(3)
      real, intent(out):: dist2
      real:: dx,dy,dz,dx1,dy1,dz1
      
      dx=abs(xyz1(1)-xyz2(1))
      dx1=box(1)-abs(xyz1(1)-xyz2(1))
      if(dx1.lt.dx) dx=dx1
      dy=abs(xyz1(2)-xyz2(2))
      dy1=box(2)-abs(xyz1(2)-xyz2(2))
      if(dy1.lt.dy) dy=dy1
      dz=abs(xyz1(3)-xyz2(3))
      dz1=box(3)-abs(xyz1(3)-xyz2(3))
      if(dz1.lt.dz) dz=dz1

      dist2=dx**2+dy**2+dz**2

      end subroutine sqdist
