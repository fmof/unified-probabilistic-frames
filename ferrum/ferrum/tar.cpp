#include "ferrum/logging.hpp"
#include "ferrum/tar.hpp"

namespace ferrum {
  CompressedTar::CompressedTar() :
    fname_(), opened_(false) {
  }
  CompressedTar::CompressedTar(const std::string& name, const CompressedTarOperation& op) :
    fname_(name), op_(op), opened_(false) {
    open(false);
  }
  CompressedTar& CompressedTar::read() {
    bool gto = false;
    switch(op_) {
    case CompressedTarOperation::UNSPECIFIED:
    case CompressedTarOperation::READ:
      op_ = CompressedTarOperation::READ;
      gto = true;
      break;
    case CompressedTarOperation::WRITE:
      WARN << "Changing operation of existing CompressedTar object from WRITE to READ";
      if(opened_) {
	close();
      }
      gto = true;
      break;
    default:
      ERROR << "Unknown operation " << op_;
    }
    if(gto) {
      open_read();
    }
    return *this;
  }
  CompressedTar& CompressedTar::write() {
    bool gto = false;
    switch(op_) {
    case CompressedTarOperation::UNSPECIFIED:
    case CompressedTarOperation::WRITE:
      op_ = CompressedTarOperation::WRITE;
      gto = true;
      break;
    case CompressedTarOperation::READ:
      WARN << "Changing operation of existing CompressedTar object from READ to WRITE";
      if(opened_) {
	close();
      }
      gto = true;
      break;
    default:
      ERROR << "Unknown operation " << op_;
    }
    if(gto) {
      open_write();
    }
    return *this;
  }
  void CompressedTar::open_read() {
    arc_ = archive_read_new();
    opened_ = true;
    archive_read_support_filter_all(arc_);
    archive_read_support_format_all(arc_);
  }
  void CompressedTar::open_write() {
    arc_ = archive_write_new();
    opened_ = true;
    archive_write_add_filter_gzip(arc_);
    archive_write_set_format_pax_restricted(arc_);
    archive_write_open_filename(arc_, fname_.c_str());
  }
  void CompressedTar::open(bool display_error) {
    switch(op_) {
    case CompressedTarOperation::READ:
      open_read();
      break;
    case CompressedTarOperation::WRITE:
      open_write();
      break;
    default:
      if(display_error) {
	ERROR << "Operation of " << op_ << " is not set properly for file " << fname_;
      }
      break;
    }
  }
  CompressedTar::~CompressedTar() {
    if(opened_) {
      close();
    }
  }
  void CompressedTar::close() {
    int r = -1;
    switch(op_) {
    case CompressedTarOperation::READ:
      r = archive_read_free(arc_);
      break;
    case CompressedTarOperation::WRITE:
      r = archive_write_close(arc_);
      archive_write_free(arc_);
      break;
    default:
      break;
    }
    if (r != ARCHIVE_OK) {
      WARN << "There was an error closing a compressed tar archive " << fname_;
    }
    opened_ = false;
  }

  int CompressedTar::write_data(const std::string& entry_name, const void* buff, size_t amount_to_write) {
    archive_entry* ae = archive_entry_new();
    archive_entry_set_pathname(ae, entry_name.c_str());
    archive_entry_set_size(ae, amount_to_write);
    archive_entry_set_filetype(ae, AE_IFREG);
    archive_entry_set_perm(ae, 0644);
    int r = archive_write_header(arc_, ae);
    CHECK_LIBARCHIVE_RESULT(r);
    archive_write_data(arc_, buff, amount_to_write);
    r |= archive_write_finish_entry(arc_);
    CHECK_LIBARCHIVE_RESULT(r);
    archive_entry_free(ae);
    return r;
  }

  CTIterator::CTIterator() : entry_(NULL), arc_(NULL) {
  }
  CTIterator::CTIterator(struct archive* a) : entry_(NULL), arc_(a) {
  }

  CTIterator&
  CTIterator::operator++() {
    int r = archive_read_next_header(arc_, &entry_);
    if( r != ARCHIVE_OK ) {
      entry_ = NULL;
    }
    return *this;
  }
  bool operator==(const CTIterator& i1, const CTIterator& i2) {
    return i1.entry_ == i2.entry_;
  }
  bool operator!=(const CTIterator& i1, const CTIterator& i2) {
    return i1.entry_ != i2.entry_;
  }
  CompressedTar::iterator
  CompressedTar::begin() {
    if(opened_) {
      close();
      open(true);
    }
    if(op_ != CompressedTarOperation::READ) {
      ERROR << "Can only begin an iterator for READ objects";
      return CTIterator();
    }
    int res = archive_read_open_filename(arc_, fname_.c_str(), 10240);
    if(res != ARCHIVE_OK) {
      throw res;
    }
    CTIterator cit(arc_);
    return ++cit;
  }
  CompressedTar::iterator
  CompressedTar::end() {
    return CTIterator();
  }

  CompressedTar::pointer
  CTIterator::operator->() {
    return entry_;
  }
  CompressedTar::pointer
  CTIterator::operator()() {
    return entry_;
  }
  CompressedTar::reference
  CTIterator::operator*() {
    return entry_;
  }

  int
  CompressedTar::operator()(const void** buff, size_t* amount_read, off_t* offset) {
    int r = archive_read_data_block(arc_, buff, amount_read, offset);
    if (r == ARCHIVE_EOF)
      return (ARCHIVE_OK);
    return (r);
  }
  int
  CompressedTar::operator()(const void** buff, size_t* size) {
    off_t o;
    return this->operator()(buff, size, &o);
  }
  void
  CompressedTar::read(const void** buff, size_t* size) {
    off_t o;
    do {
      if(this->operator()(buff, size, &o)) {
      }
    } while(0);
  }

  int
  CompressedTar::operator()(void* buff, size_t amount_to_read) {
    int r = archive_read_data(arc_, buff, amount_to_read);
    if (r == ARCHIVE_EOF)
      return (ARCHIVE_OK);
    return (r);
  }
  void
  CompressedTar::read(void* buff, size_t size) {
    do {
      if(this->operator()(buff, size)) {
      }
    } while(0);
  }
  int
  CompressedTar::operator()(const void** buff) {
    size_t s;
    off_t o;
    return this->operator()(buff, &s, &o);
  }
  std::string CompressedTar::name(CompressedTar::value_type entry) {
    return archive_entry_pathname(entry);
  }
  size_t CompressedTar::size(CompressedTar::value_type entry) {
    return archive_entry_size(entry);
  }
}
