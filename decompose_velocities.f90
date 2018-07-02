module decompose_velocities
implicit none

contains
subroutine decompose_two_velocities(coords,velocities,&
           velocity_translation,velocity_vibration,velocity_rotation)

implicit none
real,dimension(6),intent(in) :: coords,velocities
real,dimension(3),intent(out) :: velocity_translation,velocity_vibration,velocity_rotation
real,dimension(3) :: parallel_vector,orthogonal_vector
real,dimension(3) :: velocity_translation1,velocity_translation2,velocity_translation3
real,dimension(3) :: velocity_rotation1,velocity_rotation2
real,dimension(6) :: velocities_parallel
real,dimension(6) :: velocities_orthogonal
real,dimension(6) :: velocities_ortho2
real,dimension(6) :: velocities_ortho1

parallel_vector = coords(1:3) - coords(4:6)
parallel_vector = parallel_vector / sqrt(parallel_vector(1)**2 +&
					 parallel_vector(2)**2 +&
					 parallel_vector(3)**2 )

velocities_parallel(1:3) = dot_product(velocities(1:3),&
				      parallel_vector) * parallel_vector
velocities_parallel(4:6) = dot_product(velocities(4:6),&
				       parallel_vector) * parallel_vector

velocities_orthogonal = velocities - velocities_parallel

orthogonal_vector = velocities_orthogonal(1:3) + velocities_orthogonal(4:6)
orthogonal_vector = orthogonal_vector / sqrt(orthogonal_vector(1)**2 +&
					     orthogonal_vector(2)**2 +&
					     orthogonal_vector(3)**2 )

velocities_ortho1(1:3) = dot_product(velocities_orthogonal(1:3),&
				     orthogonal_vector) * orthogonal_vector
velocities_ortho1(4:6) = dot_product(velocities_orthogonal(4:6),&
				     orthogonal_vector) * orthogonal_vector

velocities_ortho2 = velocities_orthogonal - velocities_ortho1

velocity_translation1 = 0.5 * (velocities_parallel(1:3) + velocities_parallel(4:6))
velocity_vibration = velocities_parallel(1:3) - velocity_translation1(1:3)

velocity_translation2 = 0.5 * (velocities_ortho1(1:3) + velocities_ortho1(4:6))
velocity_rotation1 = velocities_ortho1(1:3) - velocity_translation2(1:3)
velocity_translation3 = 0.5 * (velocities_ortho2(1:3) + velocities_ortho2(4:6))
velocity_rotation2 = velocities_ortho2(1:3) - velocity_translation3(1:3)

velocity_translation = velocity_translation1 + &
		       velocity_translation2 + &
		       velocity_translation3
velocity_rotation = velocity_rotation1 + &
		    velocity_rotation2

end subroutine decompose_two_velocities


subroutine decompose_three_velocities(velocities,rel_velocities,velocity_center_of_mass)
implicit none
real, dimension(9),intent(in) :: velocities
real, dimension(9),intent(out) :: rel_velocities
real, dimension(3),intent(out) :: velocity_center_of_mass

velocity_center_of_mass = (velocities(1:3) + velocities(4:6) +&
			   velocities(7:9)) / 3
rel_velocities(1:3) = velocities(1:3) - velocity_center_of_mass
rel_velocities(4:6) = velocities(4:6) - velocity_center_of_mass
rel_velocities(7:9) = velocities(7:9) - velocity_center_of_mass

end subroutine decompose_three_velocities





end module decompose_velocities
