require "formula"

# This formula was created by:
# (1) Checking out the most recent formula for *thrift 0.9.1* (as
#     opposed to the most recent formula for *thrift*, since
#     thrift is currently at 0.9.2):
#       cd /usr/local/Library/Formula
#       git log thrift.rb
#       git checkout a61e5573f1e2bf06698038802b653f4a005e0743
# (2) Changing the filename to 'thrift-0.9.1.rb' and the class
#     name to Thrift091
# (3) Commenting out the "bottle" section, which provides binary
#     versions of Thrift
#
# 


class Thrift091 < Formula
  homepage "http://thrift.apache.org"

  stable do
    url "http://archive.apache.org/dist/thrift/0.9.1/thrift-0.9.1.tar.gz"
    sha1 "dc54a54f8dc706ffddcd3e8c6cd5301c931af1cc"

    # These patches are 0.9.1-specific and can go away once a newer version is released
    [
      # Apply THRIFT-2201 fix from master to 0.9.1 branch (required for clang to compile with C++11 support)
      %w{836d95f9f00be73c6936d407977796181d1a506c 4bc8c19c51f3d9f30799251a810dd1ca63c4bf1e},
      # Apply THRIFT-667
      %w{12c09f44cb291b1ecc4074cb3a55775b375fa8b2 3c449146e13727a9dbeb66e9826968b85a6d8869},
      # Apply THRIFT-1755
      %w{8cd3efe50a42975375e8ff3bc03306d9e4174314 5b45d692b2e6a4480d088485d4f011b9cd7fff52},
      # Apply THRIFT-2045
      %w{282e440c6de219b7b8f32b01cc7eb599f534f33f a5ab25621be48a803de7f0893ada1cdadaa1bc11},
      %w{9f9cd10e813ef574dd5578d78ca26a9088383d3a 65db24f59989eab8b3b0a0f1393a52faa3cba7ee},
      %w{e957675647d3d7caafe842aa85cbd987e91b21f9 478a5ad6b08a268460793112361726976fee8157},
      # Apply THRIFT-2229 fix from master to 0.9.1 branch
      %w{5f2d34e5ab33651059a085525b3adbab6a877e6f ff01eeca821a1d97485c9075c4d65f0c5eef4931},
    ].each do |name, sha|
      patch do
        url "https://git-wip-us.apache.org/repos/asf?p=thrift.git;a=patch;h=#{name}"
        sha1 sha
      end
    end
  end

#  bottle do
#    cellar :any
#    sha1 "07614d7e556b72d53e990de3966b67f8cbea88d6" => :yosemite
#    sha1 "983a86c23cc80f40f67c3897dba412e2eb3c0d97" => :mavericks
#    sha1 "57e63acf576ff07b549b9f84034d95161cb2c71c" => :mountain_lion
#  end

  head do
    url "https://git-wip-us.apache.org/repos/asf/thrift.git"

    depends_on :autoconf
    depends_on :automake
    depends_on :libtool
    depends_on "pkg-config" => :build
    depends_on "bison" => :build
  end

  option "with-haskell", "Install Haskell binding"
  option "with-erlang", "Install Erlang binding"
  option "with-java", "Install Java binding"
  option "with-perl", "Install Perl binding"
  option "with-php", "Install PHP binding"

  depends_on "boost"
  depends_on "openssl"
  depends_on :python => :optional

  def install
    system "./bootstrap.sh" unless build.stable?

    exclusions = ["--without-ruby", "--without-tests", "--without-php_extension"]

    exclusions << "--without-python" if build.without? "python"
    exclusions << "--without-haskell" if build.without? "haskell"
    exclusions << "--without-java" if build.without? "java"
    exclusions << "--without-perl" if build.without? "perl"
    exclusions << "--without-php" if build.without? "php"
    exclusions << "--without-erlang" if build.without? "erlang"

    ENV.cxx11 if MacOS.version >= :mavericks && ENV.compiler == :clang

    # Don't install extensions to /usr:
    ENV["PY_PREFIX"] = prefix
    ENV["PHP_PREFIX"] = prefix

    system "./configure", "--disable-debug",
                          "--prefix=#{prefix}",
                          "--libdir=#{lib}",
                          *exclusions
    ENV.j1
    system "make"
    system "make", "install"
  end

  def caveats
    <<-EOS.undent
    To install Ruby binding:
      gem install thrift

    To install PHP extension for e.g. PHP 5.5:
      brew install homebrew/php/php55-thrift
    EOS
  end
end
