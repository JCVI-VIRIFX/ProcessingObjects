package ProcessingObjects::SafeIO;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use Cwd 'abs_path';
use POSIX;
use Carp;

use constant DEFAULT_DIR_PERMISSIONS  => 0775;
use constant DEFAULT_FILE_PERMISSIONS => 0644;
use constant MAX_RETRY                => 3;
use constant DELAY_BETWEEN_TRIES      => 1;

our ($RETRY_NUM,$DELAY);

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT);
    $VERSION = sprintf "%d", q$Revision: 313$ =~ /(\d+)/g;
    @ISA = qw(Exporter);
    @EXPORT = qw(mkdir_safe chdir_safe rmdir_safe unlink_safe symlink_safe
                 copy_safe move_safe system_safe chmod_safe mk_tree_safe);
    $RETRY_NUM = MAX_RETRY;
    $DELAY = DELAY_BETWEEN_TRIES;
}

sub retry($@) {
    my ($func,@opts) = @_;
    my ($success,$cnt) = (0,0);
    while (!$success && $cnt < $RETRY_NUM) {
        if ($func->( @opts )) {
            $success = 1;
        } else {
            ++$cnt;
            sleep($DELAY);
        }
    }
    croak $! . ": $opts[0]" unless $success;
}

=head2 chdir_safe()

chdir_safe($path);

It tries to change working directory to the provided path. If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

=cut

sub chdir_safe($) {
    my ($path) = @_;
    retry(sub {chdir $_[0]}, $path);
}

=head2 mkdir_safe()

mkdir_safe($path);
mkdir_safe($path, $permission_mask);


It tries to create a new directory. Optionally, it is possible to provide the permission mask.
If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

Note: if you need to create more levels of directory at once or to set permission more open than the user's default mask, use mk_tree_safe() instead.

=cut

sub mkdir_safe($;$) {
    my ($file,$mode) = @_;
    retry(sub {defined $_[1] ? mkdir $_[0],$_[1] : mkdir $_[0]}, ($file,$mode));
}

=head2 rmdir_safe()

rmdir_safe($path);

It tries to remove the last directory of the provided path. If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

=cut

sub rmdir_safe($) {
    my ($file) = @_;
    retry(sub {rmdir $_[0]}, $file);
}

=head2 unlink_safe()

unlink_safe(@file_list);

It tries to remove all the files in the list. If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

=cut

sub unlink_safe(@) {
    my @fileList = @_;
    retry(sub {unlink @_}, @fileList);
}

=head2 symlink_safe()

symlink_safe($src, $dest);

It tries to create a symbolic link from the source to the destination. If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

=cut

sub symlink_safe($$) {
    my ($src,$dst) = @_;
    retry(sub {symlink $_[0],$_[1]}, $src,$dst);
}

=head2 copy_safe()

copy_safe($src, $dest);

It tries to copy a file from the source ($src) to the destination ($dest). If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

=cut

sub copy_safe($$) {
    my ($src,$dst) = @_;
    retry(sub {copy @_}, ($src,$dst));
}

=head2 move_safe()

move_safe($src, $dest);

It tries to move a file from the source ($src) to the destination ($dest). If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a fatal exception.

=cut

sub move_safe($$) {
    my ($src,$dst) = @_;
    retry(sub {move @_}, ($src,$dst));
}

=head2 chmod_safe()

chmod_safe($path, $permissions);

E.g. chmod_safe('/usr/local/scratch/tmp.txt', 0664);

It tries to change permissions of the given file or directory to the supplied values. If unsuccessful after a maximum number of tries specified by the constant MAX_RETRY, it will raise a warning, and return undef. No dieing on the spot.

=cut

sub chmod_safe($$) {
    my ($elem,$perm) = @_;
    
    unless (defined($perm)) {
        $perm = -f $elem ? DEFAULT_FILE_PERMISSIONS : DEFAULT_DIR_PERMISSIONS;
    } 
    if (eval{retry(sub {chmod, @_}, ($elem, $perm))}) {
        return(1);
    }
    else {
        carp "$@ - Possibly you aren't the owner of the file directory you're tyring to change permissions.";
        return(undef);
    }
}

=head2 system_safe()

system_safe($command_string);

It tries to execute the supplied command line. If unsuccessful at the first trial, it will raise a fatal exception.

=cut

sub system_safe($) {
    my ($cmd) = @_;
    system $cmd;
    my $er = $?;
    if (WIFEXITED($er)) {
        if (WEXITSTATUS($er) == 0) {
            # success
        } else {
            croak "system call exited with value: ".WEXITSTATUS($er);
        }
    } elsif (WIFSIGNALED($er)) { # process exited due to signal
        croak "system call exited with signal: ".WTERMSIG($er);
    } else {
        croak "system call didn't exit, and didn't get signal. You figure out what happend.";
    }
}

=item   mk_tree_safe()

my $success = mk_tree_safe($path);
my $success = mk_tree_safe($path, $permissions);

It recursively creates directories with either the default permissions of constant DEFAULT_DIR_PERMISSIONS or with the given numeric permission.
It returns the number of directories created.

E.g. mk_tree_safe('/usr/local/scratch/TMP_DIR', 0777) || die "Impossible to create the temprary directory\n\n";

=cut

sub mk_tree_safe($;$) {
    my ($dir, $perm) = @_;
    my $created = 0;
    $perm = DEFAULT_DIR_PERMISSIONS unless defined($perm);
    my @creandae = ();
    my $upper_dir = dirname($dir);
    
    while ($upper_dir ne $dir) { # Till we have elements in the path...
        last if -e $dir;
        push(@creandae, $dir);
        $dir = $upper_dir;
        $upper_dir = dirname($dir);
    }
    
    for (my $n = $#creandae; $n >= 0; --$n) {
        my $new_dir = $creandae[$n];
        mkdir_safe($new_dir);
        chmod_safe($new_dir, $perm);
        ++$created;
    }
    return($created);
}

1;
