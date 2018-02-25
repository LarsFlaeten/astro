#ifndef _ASTRO_EXCEPTIONS_H_
#define _ASTRO_EXCEPTIONS_H_

#include <exception>

namespace astro
{

class AstroException : public std::exception
{
public:
    explicit AstroException(const std::string& _short_msg, const std::string& _long_msg) : short_msg(_short_msg), long_msg(_long_msg)
    {

    }

    explicit AstroException(const char* _short_msg, const char* _long_msg) : short_msg(_short_msg), long_msg(_long_msg)
    {

    } 


	explicit AstroException(const std::string& _short_msg)
      : short_msg(_short_msg), long_msg("")
    {

    }

    explicit AstroException(const char* _short_msg)
      : short_msg(_short_msg), long_msg("")
    {

    }


	virtual const char* what() const noexcept
    {
        return short_msg.c_str();
    }

    const std::string& getShortMessage() const
    {
        return short_msg;
    }

    const std::string& getLongMessage() const
    {
        return long_msg;
    }
private:
	std::string short_msg;
    std::string long_msg;
	

};


class SpiceException : public std::exception
{
public:
    explicit SpiceException(const std::string& _short_msg, const std::string& _long_msg, const std::string& _short_msg_expl) : short_msg(_short_msg), long_msg(_long_msg), short_msg_expl(_short_msg_expl)
    {

    }
    
    explicit SpiceException(const std::string& _short_msg)
      : short_msg(_short_msg), long_msg(""), short_msg_expl("")
    {
    
    }

    virtual const char* what() const noexcept
    {
        return short_msg.c_str();
    }

    const std::string& getShortMessage() const
    {
        return short_msg;
    }
    
    const std::string& getLongMessage() const
    {
        return long_msg;
    }
    
    const std::string& getShortMessageExplanation() const
    {
        return short_msg_expl;
    }
 
private:
   // Prevent overload of this one
    explicit SpiceException(const char* what_arg)
    {}
    
    std::string short_msg;
    std::string long_msg;
    std::string short_msg_expl; 
    
};

}
#endif
