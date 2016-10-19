#ifndef FERRUM_LIBNAR_TAR_H_
#define FERRUM_LIBNAR_TAR_H_

#include <archive.h>
#include <archive_entry.h>
#include <string>

namespace ferrum {
  class CTIterator {
  public:
    CTIterator();
    CTIterator(struct archive* );
    typedef struct archive_entry* reference;
    typedef struct archive_entry* pointer;
    pointer entry_;
    pointer operator->();
    reference operator*();
    CTIterator& operator++(); // prefix
    friend bool operator==(const CTIterator&, const CTIterator&);
    friend bool operator!=(const CTIterator&, const CTIterator&);
    pointer operator()();
  private:
    struct archive* arc_; // this is shared
  };

  typedef struct archive_entry* ARC;

  enum CompressedTarOperation {
    UNSPECIFIED = 0,
    READ,
    WRITE
  };

#define CHECK_LIBARCHIVE_RESULT(r) do { \
  int ____r____ = (r);			\
  if( (____r____) != ARCHIVE_OK) {	\
  ERROR << "CompressedTar encountered an error in " << __func__ << " with file " << fname_; \
  }					\
  } while(0)

  class CompressedTar {
  public:
    typedef CTIterator iterator;
    typedef struct archive_entry* reference;
    typedef struct archive_entry* value_type;
    typedef struct archive_entry* pointer;
    CompressedTar(const std::string& name, const CompressedTarOperation& cto = CompressedTarOperation::UNSPECIFIED);
    CompressedTar();
    CompressedTar(const CompressedTar&) = delete;
    ~CompressedTar();
    CompressedTar& read();
    CompressedTar& write();
    
    iterator begin();
    iterator end();
    /**
     * do a zero-copy read into a void* buffer
     */
    int operator()(const void** buff, size_t* amount_read, off_t* offset);
    int operator()(const void** buff, size_t* amount_read);
    void read(const void** buff, size_t* amount_read);
    //////////////////////////
    int operator()(void* buff, size_t amount_to_read);
    void read(void* buff, size_t amount_to_read);

    /**
     * do a zero-copy read into a void* buffer
     */
    int operator()(const void** buff);
    std::string name(value_type entry);
    size_t size(value_type entry);
    ///////////////////////////
    int write_data(const std::string& name, const void*, size_t amount_to_write);
  private:
    std::string fname_;
    struct archive* arc_;
    //pointer entry_;
    CompressedTarOperation op_;
    bool opened_;
    void close();
    void open_read();
    void open_write();
    void open(bool disp_err);
  };
}

#endif
