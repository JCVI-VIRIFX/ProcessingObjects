package SafeIO;

use strict;
use warnings;
use File::Copy;
use POSIX;
use Carp;

our ($RETRY_NUM,$DELAY);

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT);
    $VERSION = sprintf "%d.%03d", q$Revision$ =~ /(\d+)/g;
    @ISA = qw(Exporter);
    @EXPORT = qw(mkdir_safe chdir_safe rmdir_safe unlink_safe symlink_safe
                 copy_safe move_safe system_safe);
    $RETRY_NUM = 3;
    $DELAY = 1;
}

sub retry($@) {
    my ($func,@opts) = @_;
    my ($success,$cnt) = (0,0);
    while (!$success && $cnt < $RETRY_NUM) {
        if ($func->( @opts )) {
            $success = 1;
        } else {
            $cnt++;
            sleep $DELAY;
        }
    }
    croak $! . ": $opts[0]" unless $success;
}
sub chdir_safe($) {
    my ($path) = @_;
    retry(sub {chdir $_[0]}, $path);
}
sub mkdir_safe($;$) {
    my ($file,$mode) = @_;
    retry(sub {defined $_[1] ? mkdir $_[0],$_[1] : mkdir $_[0]}, ($file,$mode));
}
sub rmdir_safe($) {
    my ($file) = @_;
    retry(sub {rmdir $_[0]}, $file);
}
sub unlink_safe(@) {
    my @fileList = @_;
    retry(sub {unlink @_}, @fileList);
}
sub symlink_safe($$) {
    my ($src,$dst) = @_;
    retry(sub {symlink $_[0],$_[1]}, $src,$dst);
}
sub copy_safe($$) {
    my ($src,$dst) = @_;
    retry(sub {copy @_}, ($src,$dst));
}
sub move_safe($$) {
    my ($src,$dst) = @_;
    retry(sub {move @_}, ($src,$dst));
}

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

1;
