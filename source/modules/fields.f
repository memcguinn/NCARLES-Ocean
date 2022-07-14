      module fields
        real, allocatable ::
     +      u(:,:,:), v(:,:,:), w(:,:,:), t(:,:,:,:), e(:,:,:),
     +      r1(:,:,:), r2(:,:,:), r3(:,:,:), r4(:,:,:,:), r5(:,:,:)
        real, allocatable ::
     +      ux(:,:,:), uy(:,:,:), vx(:,:,:), vy(:,:,:),
     +      wx(:,:,:), wy(:,:,:),
     +      p(:,:,:), ptop(:,:,:), vis_m(:,:,:), vis_s(:,:,:),
     +      vis_sv(:,:,:)
        real, allocatable ::
     +      ubc(:,:,:), vbc(:,:,:), wbc(:,:,:), tbc(:,:,:,:),
     +      ebc(:,:,:), pbc(:,:,:), pbc2(:,:,:)
      end module fields
