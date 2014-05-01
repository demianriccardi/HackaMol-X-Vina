  package HackaMol::X::Vina;
  #ABSTRACT: HackaMol extension for running Autodock Vina  
  use Moose;
  use MooseX::StrictConstructor;
  use Moose::Util::TypeConstraints;
  use Math::Vector::Real;
  use namespace::autoclean;
  use Carp;

  with qw(HackaMol::X::ExtensionRole);

  has 'receptor'  => (is => 'ro', isa => 'Str', predicate => 'has_receptor');
  has 'ligand'    => (is => 'ro', isa => 'Str', predicate => 'has_ligand');

  has $_ => (
      is => 'rw', isa => 'Num', predicate => "has_$_",
  ) foreach qw(center_x center_y center_z size_x size_y size_z);

  has $_ => (
    is => 'ro', isa => 'Int', predicate => "has_$_",
  ) foreach qw(energy_range exhaustiveness seed cpu num_modes);

  has 'center' => (
      is => 'rw', isa => 'Math::Vector::Real', predicate => "has_center",
      trigger => \&_set_center,
  );

  has 'size' => (
      is => 'rw', isa => 'Math::Vector::Real', predicate => "has_size",
      trigger => \&_set_size,
  );

  sub BUILD {
    my $self = shift;

    if ( $self->has_scratch ) {
        $self->scratch->mkpath unless ( $self->scratch->exists );
    }

#  setting center via MVR overrides whatever is set via center_x etc.
#    if ($self->has_center){
#      my ($x,$y,$z) = @{$self->center};
#      $self->center_x($x);  
#      $self->center_y($y);  
#      $self->center_z($z);  
#    }

#    if ($self->has_size){
#      my ($x,$y,$z) = @{$self->size};
#      $self->size_x($x); 
#      $self->size_y($y);
#      $self->size_z($z);
#    }
    
    unless ( $self->has_command ) {
        return unless ( $self->has_exe );
        my $cmd = $self->build_command;
        $self->command($cmd);
    }
    return;
  }

  sub _set_center {
    my ($self,$center,$old_center) = @_;
    $self->center_x($center->[0]);
    $self->center_y($center->[1]);
    $self->center_z($center->[2]);
  }

  sub _set_size {
    my ($self,$size,$old_size) = @_;
    $self->size_x($size->[0]);
    $self->size_y($size->[1]);
    $self->size_z($size->[2]);
  }

  #required methods
  sub build_command {
    my $self = shift;
    my $cmd;
    $cmd  = $self->exe;
    $cmd .= " --config " . $self->in_fn->stringify  if $self->has_in_fn;
    # we always capture output 
    return $cmd;
  }

  sub _build_map_in{
    return sub { return ( shift->write_input ) };
  }

  sub _build_map_out{
    my $sub_cr = sub {     
                      my $self = shift; 
                      my $qr = qr/^\s+\d+\s+(-*\d+\.\d)/;
                      my ($stdout,$sterr) = $self->capture_sys_command; 
                      my @be = map { m/$qr/; $1 }
                               grep{ m/$qr/ } 
                               split ("\n",$stdout);  
                      return (@be);
                     };
    return $sub_cr;
  }

  sub write_input {
    my $self  = shift;

    unless ($self->has_in_fn) {
      croak "no vina in_fn for writing input";
    }

    my $input ;
    $input   .= sprintf("%-15s = %-55s\n",'out', $self->out_fn->stringify) if $self->has_out_fn;
    $input   .= sprintf("%-15s = %-55s\n",'log', $self->log_fn->stringify) if $self->has_log_fn;
    foreach my $cond (qw(receptor ligand cpu num_modes energy_range exhaustiveness seed)) {
      my $condition = "has_$cond";
      $input .= sprintf("%-15s = %-55s\n",$cond , $self->$cond) if $self->$condition;
    }
    foreach my $metric (qw(center_x center_y center_z size_x size_y size_z)) {
      $input .= sprintf("%-15s = %-55s\n",$metric , $self->$metric);
    }
    $self->in_fn->spew($input);
    return ($input); 
  }

  __PACKAGE__->meta->make_immutable;

  1;

__END__

=head1 SYNOPSIS

  

